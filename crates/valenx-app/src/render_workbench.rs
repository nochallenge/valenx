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

impl RenderWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]): the four numeric render
    /// parameters plus the stainless-finish toggle.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "resolution",
            "samples / px",
            "max bounces",
            "exposure",
            "polished stainless finish (rocket)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. The three count fields read `AgentValue::as_i64`
    /// (the render clamps the value at solve time, so any non-negative integer
    /// is accepted here); `exposure` reads `AgentValue::as_f64` (narrowed to
    /// `f32`); the finish toggle reads `AgentValue::as_bool`. Fail-loud: an
    /// unknown caption, wrong type, or negative count returns `Err` — never a
    /// panic, no field written on error.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "resolution" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("resolution must be >= 0, got {n}"));
                }
                self.resolution = n as u32;
            }
            "samples / px" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("samples / px must be >= 0, got {n}"));
                }
                self.spp = n as u32;
            }
            "max bounces" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("max bounces must be >= 0, got {n}"));
                }
                self.max_depth = n as u32;
            }
            "exposure" => self.exposure = value.as_f64()? as f32,
            "polished stainless finish (rocket)" => self.stainless = value.as_bool()?,
            other => return Err(format!("unknown Render control: {other:?}")),
        }
        Ok(())
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
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_render_workbench",
        "Path-Traced Render",
        |app, ui| {
            ui.label(
                egui::RichText::new("global illumination · valenx-pathtrace")
                    .weak()
                    .small(),
            );
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
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as its
                    // accessibility / UI-Automation Name (egui clears a DragValue's
                    // own Name, leaving it anonymous to a screen reader / AI driver
                    // otherwise).
                    let res = ui.label("resolution");
                    ui.add(
                        egui::DragValue::new(&mut s.resolution)
                            .speed(4.0)
                            .range(48..=512),
                    )
                    .labelled_by(res.id);
                    ui.end_row();
                    let spp = ui.label("samples / px");
                    ui.add(egui::DragValue::new(&mut s.spp).speed(1.0).range(1..=128))
                        .labelled_by(spp.id);
                    ui.end_row();
                    let mb = ui.label("max bounces");
                    ui.add(
                        egui::DragValue::new(&mut s.max_depth)
                            .speed(0.2)
                            .range(1..=16),
                    )
                    .labelled_by(mb.id);
                    ui.end_row();
                    let exp = ui.label("exposure");
                    ui.add(
                        egui::DragValue::new(&mut s.exposure)
                            .speed(0.05)
                            .range(0.1..=4.0),
                    )
                    .labelled_by(exp.id);
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
        },
    );
    if close {
        app.show_render_workbench = false;
    }

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

/// Resolution (square) of the agent-bridge `render` product's path-trace.
/// Small on purpose — the registry builder is **synchronous** (it runs inline
/// when the `show_3d`/registry reducer resolves the kind, with no background
/// worker), so the render must finish in a few milliseconds. 96² with a
/// handful of samples is enough to read as a recognisable Cornell box in the
/// tile.
const PRODUCT_RENDER_RES: u32 = 96;
/// Samples-per-pixel for the agent-bridge `render` product — low so the
/// synchronous build stays fast (a noisy-but-legible preview, like the `cfd`
/// card's bounded solve).
const PRODUCT_RENDER_SPP: u32 = 8;
/// Path-trace bounce depth for the agent-bridge `render` product.
const PRODUCT_RENDER_DEPTH: u32 = 4;
/// Exposure for the agent-bridge `render` product's tone-map.
const PRODUCT_RENDER_EXPOSURE: f32 = 1.0;

/// Build the agent-bridge **`render` product** — a small synchronous Cornell-box
/// path-trace whose tone-mapped RGB8 framebuffer becomes the tile's
/// [`egui::ColorImage`] (the IMAGE product the pane renders scaled-to-fit). The
/// matching `image: Some(..)` path in [`crate::dock_layout`]'s
/// `render_workspace_body` uploads it to a texture once and draws it.
///
/// Unlike the interactive panel (which renders on a worker thread and uploads a
/// texture there), this builder is **pure and app-state-free** so the registry
/// reducer can call it inline; the resolution / samples are kept small
/// ([`PRODUCT_RENDER_RES`] / [`PRODUCT_RENDER_SPP`]) so that inline render is a
/// few-millisecond operation rather than the panel's full-quality render. The
/// `lines` carry the render parameters as readout rows. The `48..=512`
/// resolution clamp inside [`render_demo`] keeps its only error
/// (`FramebufferError::TooLarge`) unreachable at this size; if it ever did
/// fail, the product degrades to a parameter-only text card (no `image`) rather
/// than panicking on the bridge path.
pub(crate) fn render_product() -> crate::WorkspaceProduct {
    let lines = match render_demo(
        PRODUCT_RENDER_RES,
        PRODUCT_RENDER_SPP,
        PRODUCT_RENDER_DEPTH,
        PRODUCT_RENDER_EXPOSURE,
    ) {
        Ok((w, h, pixels)) => {
            // RGB8 → egui ColorImage. `from_rgb` expects exactly `w·h·3` bytes,
            // which `render_demo` guarantees; the texture upload happens lazily
            // in the tile renderer.
            let image = egui::ColorImage::from_rgb([w, h], &pixels);
            return crate::WorkspaceProduct {
                title: "Path-traced render".into(),
                lines: vec![
                    format!("Cornell box · {w}×{h} px"),
                    format!("{PRODUCT_RENDER_SPP} samples/px · {PRODUCT_RENDER_DEPTH} bounces"),
                    "valenx-pathtrace · global illumination".into(),
                ],
                mesh: None,
                vertex_colors: None,
                camera: valenx_viz::OrbitCamera::default(),
                kind2d: None,
                last_export: None,
                image: Some(image),
                image_texture: None,
                animation: None,
            };
        }
        // Unreachable at 96² (well under the cap) — but the bridge path must not
        // panic, so degrade to a text card describing the intended render.
        Err(e) => vec![
            "Path-traced render (Cornell box)".into(),
            format!("render unavailable: {e}"),
        ],
    };
    crate::WorkspaceProduct {
        title: "Path-traced render".into(),
        lines,
        mesh: None,
        vertex_colors: None,
        camera: valenx_viz::OrbitCamera::default(),
        kind2d: None,
        last_export: None,
        image: None,
        image_texture: None,
        animation: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        use crate::agent_commands::AgentValue;
        let mut s = RenderWorkbenchState::default();
        // A representative integer set lands in state.
        s.agent_set("resolution", &AgentValue::Int(256)).unwrap();
        assert_eq!(s.resolution, 256);
        // A float param narrows to f32.
        s.agent_set("exposure", &AgentValue::Float(1.5)).unwrap();
        assert_eq!(s.exposure, 1.5);
        // The checkbox is a bool toggle.
        s.agent_set(
            "polished stainless finish (rocket)",
            &AgentValue::Bool(true),
        )
        .unwrap();
        assert!(s.stainless);
        // Unknown caption -> Err.
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into the integer resolution) -> Err, untouched.
        assert!(s
            .agent_set("resolution", &AgentValue::Str("big".into()))
            .is_err());
        assert_eq!(s.resolution, 256, "rejected set leaves field untouched");
    }

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

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    /// As the panel draw, but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_render_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_render_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_render_workbench(&mut app, ctx);
        });
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The render parameters (resolution, samples/px, max bounces, exposure)
        // are SpinButtons; each must be `labelled_by` its caption (egui clears a
        // DragValue's own Name) so an AI / screen reader can find the control by
        // the caption text.
        let mut app = ValenxApp::default();
        app.show_render_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected the render numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every render DragValue must be labelled_by its caption (AI-drivable name)"
        );
        for caption in ["resolution", "max bounces"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }
}
