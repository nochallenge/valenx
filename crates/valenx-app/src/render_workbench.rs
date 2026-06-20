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

use crate::background::{BackgroundJob, JobState};
use crate::ValenxApp;

/// The path-trace worker's output: `(width, height, RGB8 pixels)` on success,
/// or a display-string error.
type RenderOutput = Result<(usize, usize, Vec<u8>), String>;

/// Persistent state for the render workbench.
pub struct RenderWorkbenchState {
    resolution: u32,
    spp: u32,
    max_depth: u32,
    exposure: f32,
    texture: Option<egui::TextureHandle>,
    status: String,
    error: Option<String>,
    /// Render the rocket in polished-stainless (vs white-painted) finish.
    stainless: bool,
    /// A render running on a worker thread. While `Some`, the form is frozen
    /// and a spinner shows; `poll_render` uploads the texture when it finishes.
    job: Option<BackgroundJob<RenderOutput>>,
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
            error: None,
            stainless: false,
            job: None,
        }
    }
}

/// Render a Cornell-box scene and return `(width, height, RGB8 pixels)`.
///
/// Returns the `valenx-pathtrace` framebuffer error as a display string
/// rather than panicking: the only failure is `FramebufferError::TooLarge`,
/// which the `48..=512` resolution clamp keeps unreachable today — but a
/// fallible call must never `expect`/panic on the user-clickable render
/// path, so the error is surfaced in the panel instead.
fn render_demo(resolution: u32, spp: u32, max_depth: u32, exposure: f32) -> RenderOutput {
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
    b.add_quad(
        vec3(-1.0, 0.0, -1.0),
        vec3(1.0, 0.0, -1.0),
        vec3(1.0, 0.0, 1.0),
        vec3(-1.0, 0.0, 1.0),
        white,
    );
    b.add_quad(
        vec3(-1.0, 2.0, -1.0),
        vec3(-1.0, 2.0, 1.0),
        vec3(1.0, 2.0, 1.0),
        vec3(1.0, 2.0, -1.0),
        white,
    );
    b.add_quad(
        vec3(-1.0, 0.0, -1.0),
        vec3(-1.0, 2.0, -1.0),
        vec3(1.0, 2.0, -1.0),
        vec3(1.0, 0.0, -1.0),
        white,
    );
    b.add_quad(
        vec3(-1.0, 0.0, -1.0),
        vec3(-1.0, 0.0, 1.0),
        vec3(-1.0, 2.0, 1.0),
        vec3(-1.0, 2.0, -1.0),
        red,
    );
    b.add_quad(
        vec3(1.0, 0.0, -1.0),
        vec3(1.0, 2.0, -1.0),
        vec3(1.0, 2.0, 1.0),
        vec3(1.0, 0.0, 1.0),
        green,
    );
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
    let ldr = render(&scene, &params)
        .map_err(|e| e.to_string())?
        .to_ldr(exposure);
    Ok((ldr.width as usize, ldr.height as usize, ldr.pixels))
}

/// Path-trace the **Valenx LV-1** as a "final picture": the procedural rocket
/// mesh in a brushed-metal finish on a ground plane, lit by a big overhead
/// light — the photoreal counterpart to the live shaded viewport.
pub(crate) fn render_rocket(
    resolution: u32,
    spp: u32,
    max_depth: u32,
    exposure: f32,
    stainless: bool,
) -> RenderOutput {
    let res = resolution.clamp(48, 512);
    // The rocket stands along +Z, so the camera's up is +Z; a 3/4 hero framing
    // of the ~33-tall vehicle.
    let camera = PtCamera::look_at(
        vec3(42.0, 38.0, 22.0),
        vec3(0.0, 0.0, 13.0),
        vec3(0.0, 0.0, 1.0),
        34f32.to_radians(),
        res,
        res,
    );
    let mut b = SceneBuilder::new(camera);
    let body = b.add_material(if stainless {
        PtMaterial::metal([0.96, 0.96, 0.97], 0.07) // polished stainless steel
    } else {
        PtMaterial::metal([0.92, 0.92, 0.93], 0.34) // white-painted launcher
    });
    let ground = b.add_material(PtMaterial::diffuse([0.46, 0.45, 0.43]));
    let key = b.add_material(PtMaterial::emissive([13.5, 12.6, 11.0]));
    let fill = b.add_material(PtMaterial::emissive([2.4, 2.4, 2.6]));

    let rocket = crate::rocket_mesh::lv1_rocket_mesh();
    b.add_mesh(&rocket, body);

    // Ground plane at the engine-bell base.
    let (g, z0) = (140.0, -2.8);
    b.add_quad(
        vec3(-g, -g, z0),
        vec3(g, -g, z0),
        vec3(g, g, z0),
        vec3(-g, g, z0),
        ground,
    );
    // Big overhead key light (emits downward, like the Cornell ceiling light).
    let (lx, ly, lz) = (50.0, 50.0, 62.0);
    b.add_quad(
        vec3(-lx, -ly, lz),
        vec3(lx, -ly, lz),
        vec3(lx, ly, lz),
        vec3(-lx, ly, lz),
        key,
    );
    // Cool sky-fill backdrop so shadowed faces aren't pure black.
    b.add_quad(
        vec3(-140.0, 95.0, -10.0),
        vec3(140.0, 95.0, -10.0),
        vec3(140.0, 95.0, 85.0),
        vec3(-140.0, 95.0, 85.0),
        fill,
    );

    let scene = b.build();
    let params = RenderParams {
        samples_per_pixel: spp.clamp(1, 128),
        max_depth: max_depth.clamp(1, 16),
        seed: 0x5eed,
        exposure,
    };
    let ldr = render(&scene, &params)
        .map_err(|e| e.to_string())?
        .to_ldr(exposure);
    Ok((ldr.width as usize, ldr.height as usize, ldr.pixels))
}

/// Path-trace the **detailed engine** (chamber, fluted cooling-channel nozzle,
/// injector dome) close-up in a metallic finish — the "what does the engine
/// look like" view.
pub(crate) fn render_engine(
    resolution: u32,
    spp: u32,
    max_depth: u32,
    exposure: f32,
) -> RenderOutput {
    let res = resolution.clamp(48, 512);
    // Frame the full engine including the tall powerhead (exit at z=0 up to
    // the turbine tops near z=12.4); a 3/4 hero angle, +Z up.
    let camera = PtCamera::look_at(
        vec3(21.0, 18.0, 10.5),
        vec3(0.0, 0.0, 6.0),
        vec3(0.0, 0.0, 1.0),
        44f32.to_radians(),
        res,
        res,
    );
    let mut b = SceneBuilder::new(camera);
    let metal = b.add_material(PtMaterial::metal([0.86, 0.86, 0.87], 0.16));
    let ground = b.add_material(PtMaterial::diffuse([0.34, 0.34, 0.36]));
    let key = b.add_material(PtMaterial::emissive([14.0, 13.2, 11.6]));
    let fill = b.add_material(PtMaterial::emissive([2.6, 2.6, 2.8]));

    let engine = crate::rocket_mesh::detailed_engine_mesh();
    b.add_mesh(&engine, metal);
    let (g, z0) = (60.0, -0.3);
    b.add_quad(
        vec3(-g, -g, z0),
        vec3(g, -g, z0),
        vec3(g, g, z0),
        vec3(-g, g, z0),
        ground,
    );
    let (lx, ly, lz) = (25.0, 25.0, 28.0);
    b.add_quad(
        vec3(-lx, -ly, lz),
        vec3(lx, -ly, lz),
        vec3(lx, ly, lz),
        vec3(-lx, ly, lz),
        key,
    );
    b.add_quad(
        vec3(-60.0, 45.0, -5.0),
        vec3(60.0, 45.0, -5.0),
        vec3(60.0, 45.0, 40.0),
        vec3(-60.0, 45.0, 40.0),
        fill,
    );

    let scene = b.build();
    let params = RenderParams {
        samples_per_pixel: spp.clamp(1, 128),
        max_depth: max_depth.clamp(1, 16),
        seed: 0x5eed,
        exposure,
    };
    let ldr = render(&scene, &params)
        .map_err(|e| e.to_string())?
        .to_ldr(exposure);
    Ok((ldr.width as usize, ldr.height as usize, ldr.pixels))
}

/// Move a finished background render into the panel: build the egui texture on
/// the UI thread (the worker only produced raw pixels), or surface the error.
fn poll_render(s: &mut RenderWorkbenchState, ctx: &egui::Context) {
    match s.job.as_mut().map(BackgroundJob::poll) {
        Some(JobState::Done(result)) => {
            s.job = None;
            match result {
                Ok((w, h, pixels)) => {
                    let mut rgba = Vec::with_capacity(w * h * 4);
                    for px in pixels.chunks_exact(3) {
                        rgba.extend_from_slice(&[px[0], px[1], px[2], 255]);
                    }
                    let color = egui::ColorImage::from_rgba_unmultiplied([w, h], &rgba);
                    let tex =
                        ctx.load_texture("pathtrace_render", color, egui::TextureOptions::LINEAR);
                    s.texture = Some(tex);
                    s.status = format!("rendered {w}×{h} @ {} spp", s.spp);
                }
                Err(e) => {
                    s.error = Some(format!("render failed: {e}"));
                }
            }
        }
        Some(JobState::Failed) => {
            s.job = None;
            s.error = Some("the render thread stopped unexpectedly".into());
        }
        Some(JobState::Pending) | None => {}
    }
}

/// Draw the render workbench (a no-op unless toggled on via View → Render).
pub fn draw_render_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_render_workbench {
        return;
    }
    poll_render(&mut app.render, ctx);
    let mut do_render = false;
    let mut do_render_rocket = false;
    let mut do_render_engine = false;
    let close = crate::workbench_chrome::workbench_shell(app, ctx, "valenx_render_workbench", "Path-Traced Render", |app, ui| {
            ui.label(egui::RichText::new("global illumination · valenx-pathtrace").weak().small());
            ui.separator();
            let s = &mut app.render;
            let running = s.job.is_some();
            if running {
                ui.ctx().request_repaint();
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("rendering…");
                });
                ui.disable();
            }
            egui::Grid::new("render_params")
                .num_columns(2)
                .show(ui, |ui| {
                    ui.label("resolution");
                    ui.add(
                        egui::DragValue::new(&mut s.resolution)
                            .speed(4.0)
                            .range(48..=512),
                    );
                    ui.end_row();
                    ui.label("samples / px");
                    ui.add(egui::DragValue::new(&mut s.spp).speed(1.0).range(1..=128));
                    ui.end_row();
                    ui.label("max bounces");
                    ui.add(
                        egui::DragValue::new(&mut s.max_depth)
                            .speed(0.2)
                            .range(1..=16),
                    );
                    ui.end_row();
                    ui.label("exposure");
                    ui.add(
                        egui::DragValue::new(&mut s.exposure)
                            .speed(0.05)
                            .range(0.1..=4.0),
                    );
                    ui.end_row();
                });
            if ui.button("▶ Render (Cornell box)").clicked() {
                do_render = true;
            }
            if ui
                .button(egui::RichText::new("🚀 Render the rocket (final picture)").strong())
                .clicked()
            {
                do_render_rocket = true;
            }
            if ui
                .button(egui::RichText::new("🔧 Render the engine (detail)").strong())
                .clicked()
            {
                do_render_engine = true;
            }
            ui.checkbox(&mut s.stainless, "polished stainless finish (rocket)");
            if !s.status.is_empty() {
                ui.label(egui::RichText::new(&s.status).small().weak());
            }
            if let Some(err) = &s.error {
                ui.colored_label(egui::Color32::from_rgb(220, 90, 90), err);
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
        }, );
    if close { app.show_render_workbench = false; }

    if do_render || do_render_rocket || do_render_engine {
        // Clear stale per-run state up front (house style — cf. cfd_workbench),
        // so a failed render can't leave the previous success line showing next
        // to the error. The last-good texture is deliberately retained.
        let s = &mut app.render;
        s.status.clear();
        s.error = None;
        let (res, spp, max_depth, exposure) = (s.resolution, s.spp, s.max_depth, s.exposure);
        let (rocket, engine, stainless) = (do_render_rocket, do_render_engine, s.stainless);
        // Render on a worker thread; poll_render uploads the texture when it
        // finishes, so the path tracer no longer freezes the UI.
        s.job = Some(BackgroundJob::spawn(move || {
            if engine {
                render_engine(res, spp, max_depth, exposure)
            } else if rocket {
                render_rocket(res, spp, max_depth, exposure, stainless)
            } else {
                render_demo(res, spp, max_depth, exposure)
            }
        }));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_a_lit_image() {
        // Small + cheap: the scene must produce a correctly-sized RGB buffer
        // with some non-black pixels (the emissive light illuminates it).
        let (w, h, pixels) =
            render_demo(64, 4, 4, 1.0).expect("64² render is well under the pixel cap");
        assert_eq!(pixels.len(), w * h * 3, "RGB8 buffer of the right size");
        assert!(
            pixels.iter().any(|&p| p > 0),
            "the scene is lit (non-black pixels)"
        );
    }

    #[test]
    fn renders_the_rocket_lit() {
        // The rocket "final picture" must produce a correctly-sized buffer
        // that is actually lit — a fair number of bright pixels, not a near-
        // black frame (which would mean the lighting/camera is wrong).
        let (w, h, pixels) =
            render_rocket(64, 8, 4, 1.0, false).expect("64² rocket render is under the pixel cap");
        assert_eq!(pixels.len(), w * h * 3, "RGB8 buffer of the right size");
        let bright = pixels.iter().filter(|&&p| p > 30).count();
        assert!(
            bright > pixels.len() / 20,
            "rocket scene should be well lit; only {bright} bright sub-pixels"
        );
    }

    /// Render the rocket at preview quality and write a PNG to TEMP — run with
    /// `cargo test -p valenx-app --lib render_workbench::tests::dump -- --ignored --nocapture`.
    #[test]
    #[ignore = "writes path-traced rocket PNGs (white + stainless) to TEMP"]
    fn dump_rocket_png() {
        for (name, stainless) in [("white", false), ("stainless", true)] {
            let (w, h, pixels) = render_rocket(460, 256, 6, 1.1, stainless).expect("rocket render");
            let path = std::env::temp_dir().join(format!("valenx_rocket_{name}.png"));
            let file = std::fs::File::create(&path).expect("create png");
            let mut enc = png::Encoder::new(std::io::BufWriter::new(file), w as u32, h as u32);
            enc.set_color(png::ColorType::Rgb);
            enc.set_depth(png::BitDepth::Eight);
            let mut writer = enc.write_header().expect("png header");
            writer.write_image_data(&pixels).expect("png data");
            writer.finish().expect("png finish");
            println!("WROTE {}", path.display());
        }
    }

    #[test]
    fn background_render_delivers_pixels() {
        // Exercise the reactive path: render_demo runs on a worker thread and
        // the result is polled — the panel no longer blocks the UI thread.
        let mut job = BackgroundJob::spawn(|| render_demo(64, 4, 4, 1.0));
        let mut out = None;
        for _ in 0..2000 {
            match job.poll() {
                JobState::Pending => std::thread::sleep(std::time::Duration::from_millis(1)),
                JobState::Done(r) => {
                    out = Some(r);
                    break;
                }
                JobState::Failed => panic!("render worker should not fail"),
            }
        }
        let (w, h, pixels) = out
            .expect("render delivered within the poll budget")
            .expect("64² render is well under the pixel cap");
        assert_eq!(pixels.len(), w * h * 3);
        assert!(pixels.iter().any(|&p| p > 0), "the scene is lit");
    }

    #[test]
    fn renders_the_engine_lit() {
        let (w, h, pixels) =
            render_engine(64, 8, 4, 1.0).expect("64² engine render is under the cap");
        assert_eq!(pixels.len(), w * h * 3);
        let bright = pixels.iter().filter(|&&p| p > 30).count();
        assert!(bright > pixels.len() / 20, "engine scene should be lit");
    }

    #[test]
    #[ignore = "writes a path-traced engine PNG to TEMP"]
    fn dump_engine_png() {
        let (w, h, pixels) = render_engine(480, 224, 6, 1.1).expect("engine render");
        let path = std::env::temp_dir().join("valenx_engine.png");
        let file = std::fs::File::create(&path).expect("create png");
        let mut enc = png::Encoder::new(std::io::BufWriter::new(file), w as u32, h as u32);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        let mut writer = enc.write_header().expect("png header");
        writer.write_image_data(&pixels).expect("png data");
        writer.finish().expect("png finish");
        println!("WROTE {}", path.display());
    }
}
