//! Path-traced render viewer — global illumination on `valenx-pathtrace`.
//!
//! A right-side panel that renders a small Cornell-box scene with the native
//! CPU path tracer and displays the result as an image (an egui texture). The
//! render is the first **image-display** UI in the app: build a scene → render
//! to an LDR pixel buffer → upload as a texture → show it. Render parameters
//! (resolution, samples-per-pixel, bounces, exposure) are live. The render
//! itself is headless-testable; the texture display is the interactive view.

use eframe::egui;

use valenx_pathtrace::{render, vec3, PtCamera, PtMaterial, RenderParams, SceneBuilder};

use crate::ValenxApp;

/// Persistent state for the render workbench.
pub struct RenderWorkbenchState {
    resolution: u32,
    spp: u32,
    max_depth: u32,
    exposure: f32,
    texture: Option<egui::TextureHandle>,
    status: String,
}

impl Default for RenderWorkbenchState {
    fn default() -> Self {
        Self {
            resolution: 192,
            spp: 12,
            max_depth: 6,
            exposure: 1.0,
            texture: None,
            status: String::new(),
        }
    }
}

/// Render a Cornell-box scene and return `(width, height, RGB8 pixels)`.
fn render_demo(resolution: u32, spp: u32, max_depth: u32, exposure: f32) -> (usize, usize, Vec<u8>) {
    let res = resolution.clamp(48, 512);
    let camera = PtCamera::look_at(
        vec3(0.0, 1.0, 3.6),
        vec3(0.0, 1.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        50f32.to_radians(),
        res,
        res,
    );
    let mut b = SceneBuilder::new(camera);
    let white = b.add_material(PtMaterial::diffuse([0.75, 0.75, 0.75]));
    let red = b.add_material(PtMaterial::diffuse([0.70, 0.15, 0.15]));
    let green = b.add_material(PtMaterial::diffuse([0.15, 0.70, 0.15]));
    let light = b.add_material(PtMaterial::emissive([8.0, 8.0, 8.0]));
    // Cornell-ish box: floor, ceiling, back wall, red left, green right.
    b.add_quad(vec3(-1.0, 0.0, -1.0), vec3(1.0, 0.0, -1.0), vec3(1.0, 0.0, 1.0), vec3(-1.0, 0.0, 1.0), white);
    b.add_quad(vec3(-1.0, 2.0, -1.0), vec3(-1.0, 2.0, 1.0), vec3(1.0, 2.0, 1.0), vec3(1.0, 2.0, -1.0), white);
    b.add_quad(vec3(-1.0, 0.0, -1.0), vec3(-1.0, 2.0, -1.0), vec3(1.0, 2.0, -1.0), vec3(1.0, 0.0, -1.0), white);
    b.add_quad(vec3(-1.0, 0.0, -1.0), vec3(-1.0, 0.0, 1.0), vec3(-1.0, 2.0, 1.0), vec3(-1.0, 2.0, -1.0), red);
    b.add_quad(vec3(1.0, 0.0, -1.0), vec3(1.0, 2.0, -1.0), vec3(1.0, 2.0, 1.0), vec3(1.0, 0.0, 1.0), green);
    // Ceiling light.
    b.add_quad(
        vec3(-0.3, 1.98, -0.3),
        vec3(0.3, 1.98, -0.3),
        vec3(0.3, 1.98, 0.3),
        vec3(-0.3, 1.98, 0.3),
        light,
    );
    let scene = b.build();
    let params = RenderParams {
        samples_per_pixel: spp.clamp(1, 128),
        max_depth: max_depth.clamp(1, 16),
        seed: 0x5eed,
        exposure,
    };
    let ldr = render(&scene, &params).expect("render succeeds").to_ldr(exposure);
    (ldr.width as usize, ldr.height as usize, ldr.pixels)
}

/// Draw the render workbench (a no-op unless toggled on via View → Render).
pub fn draw_render_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_render_workbench {
        return;
    }
    let mut do_render = false;
    egui::SidePanel::right("valenx_render_workbench")
        .resizable(true)
        .default_width(400.0)
        .width_range(340.0..=760.0)
        .show(ctx, |ui| {
            ui.heading("Path-Traced Render");
            ui.label(
                egui::RichText::new("global illumination · valenx-pathtrace")
                    .weak()
                    .small(),
            );
            ui.separator();
            let s = &mut app.render;
            egui::Grid::new("render_params").num_columns(2).show(ui, |ui| {
                ui.label("resolution");
                ui.add(egui::DragValue::new(&mut s.resolution).speed(4.0).range(48..=512));
                ui.end_row();
                ui.label("samples / px");
                ui.add(egui::DragValue::new(&mut s.spp).speed(1.0).range(1..=128));
                ui.end_row();
                ui.label("max bounces");
                ui.add(egui::DragValue::new(&mut s.max_depth).speed(0.2).range(1..=16));
                ui.end_row();
                ui.label("exposure");
                ui.add(egui::DragValue::new(&mut s.exposure).speed(0.05).range(0.1..=4.0));
                ui.end_row();
            });
            if ui.button("▶ Render (Cornell box)").clicked() {
                do_render = true;
            }
            if !s.status.is_empty() {
                ui.label(egui::RichText::new(&s.status).small().weak());
            }
            if let Some(tex) = &s.texture {
                ui.add_space(4.0);
                ui.add(
                    egui::Image::new(egui::load::SizedTexture::new(tex.id(), tex.size_vec2()))
                        .max_width(ui.available_width()),
                );
            } else {
                ui.label(
                    egui::RichText::new("Press Render — a few hundred ms for the default size.")
                        .small()
                        .weak(),
                );
            }
        });

    if do_render {
        let (w, h, pixels) =
            render_demo(app.render.resolution, app.render.spp, app.render.max_depth, app.render.exposure);
        let mut rgba = Vec::with_capacity(w * h * 4);
        for px in pixels.chunks_exact(3) {
            rgba.extend_from_slice(&[px[0], px[1], px[2], 255]);
        }
        let color = egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba);
        let tex = ctx.load_texture("pathtrace_render", color, egui::TextureOptions::LINEAR);
        app.render.texture = Some(tex);
        app.render.status = format!("rendered {w}×{h} @ {} spp", app.render.spp);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_a_lit_image() {
        // Small + cheap: the scene must produce a correctly-sized RGB buffer
        // with some non-black pixels (the emissive light illuminates it).
        let (w, h, pixels) = render_demo(64, 4, 4, 1.0);
        assert_eq!(pixels.len(), w * h * 3, "RGB8 buffer of the right size");
        assert!(pixels.iter().any(|&p| p > 0), "the scene is lit (non-black pixels)");
    }
}
