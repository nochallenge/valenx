//! The right-side **Acoustics Workbench** panel — native rectangular-room
//! reverberation + source-level analysis over `valenx-acoustics`.
//!
//! Mirrors the Heat Transfer / Gearbox workbenches: a resizable
//! [`egui::SidePanel`] gated on `crate::ValenxApp::show_acoustics_workbench`,
//! toggled from the View menu. The form sets a rectangular room (`Lx·Ly·Lz`),
//! a uniform surface absorption coefficient and the air temperature, plus a
//! point source level and a listener distance; "Analyze" reports the
//! statistical (diffuse-field) reverberation time `RT60` by both
//! [`Sabine`](valenx_acoustics::sabine_reverberation_time) and
//! [`Eyring`](valenx_acoustics::eyring_reverberation_time), the lowest axial
//! room mode, and the free-field sound-pressure level at the listener (the
//! inverse-square distance drop), and "Show 3-D room" loads a representative
//! room shell with the source inside into the central viewport.

use std::path::PathBuf;

use eframe::egui;
use nalgebra::Vector3;

use valenx_acoustics::radiation::{monopole_pressure, RHO_AIR};
use valenx_acoustics::room::RoomDimensions;
use valenx_acoustics::{
    eyring_reverberation_time, pressure_from_spl, sabine_reverberation_time, speed_of_sound, spl,
    total_absorption,
};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::types::LoadedMesh;
use crate::ValenxApp;

/// Persistent form + result state for the Acoustics Workbench.
pub struct AcousticsWorkbenchState {
    /// Room length along x `Lx` (m).
    length_x_m: f64,
    /// Room length along y `Ly` (m).
    length_y_m: f64,
    /// Room length along z (height) `Lz` (m).
    length_z_m: f64,
    /// Uniform surface absorption coefficient `ᾱ` (`0 < ᾱ < 1`).
    absorption: f64,
    /// Air temperature (deg C) — sets the speed of sound.
    temperature_c: f64,
    /// Point-source sound-pressure level at the reference distance (dB).
    source_db: f64,
    /// Reference distance at which `source_db` is quoted (m).
    ref_distance_m: f64,
    /// Listener distance from the source (m).
    listener_distance_m: f64,

    // -- Monopole radiation inputs (the pulsating-sphere source model) ----
    /// Source sphere radius `a` (m).
    source_radius_m: f64,
    /// Sphere surface (radial) velocity amplitude `U` (m/s).
    surface_velocity: f64,
    /// Radiation frequency `f` (Hz).
    frequency_hz: f64,
    /// Observer (field-point) distance `r` from the source centre (m).
    observer_distance_m: f64,
    /// Latest monopole radiation readout `(pressure_pa, spl_db)`, or `None`
    /// before the first radiation compute.
    radiation: Option<(f64, f64)>,

    /// Formatted performance readout (empty until the first analyze).
    result: String,
    /// Validation / compute error, if any.
    error: Option<String>,
    /// Deferred request to build the 3-D room solid (serviced after the
    /// panel draws).
    show_3d_request: bool,
}

impl Default for AcousticsWorkbenchState {
    fn default() -> Self {
        // A 7 x 5 x 3 m room (V = 105 m^3, S = 142 m^2) with a uniform
        // absorption coefficient of 0.15 at 20 C (c ~ 343 m/s): A = 21.3
        // sabins, Sabine RT60 ~ 0.79 s, Eyring ~ 0.73 s, lowest axial mode
        // ~ 24.5 Hz. An 85 dB source at 1 m reads ~73 dB at the 4 m listener
        // (a 12 dB inverse-square drop).
        Self {
            length_x_m: 7.0,
            length_y_m: 5.0,
            length_z_m: 3.0,
            absorption: 0.15,
            temperature_c: 20.0,
            source_db: 85.0,
            ref_distance_m: 1.0,
            listener_distance_m: 4.0,
            // Monopole radiation seed: a 0.1 m sphere pulsating at 0.01 m/s,
            // 1000 Hz, observed at 1 m in air (rho = 1.204, c = 343) — the
            // crate's validated doc example (pressure at 2 m is half that at
            // 1 m, the 1/r free-field monopole law).
            source_radius_m: 0.1,
            surface_velocity: 0.01,
            frequency_hz: 1000.0,
            observer_distance_m: 1.0,
            radiation: None,
            result: String::new(),
            error: None,
            show_3d_request: false,
        }
    }
}

/// Draw the Acoustics Workbench right-side panel. A no-op when the
/// `show_acoustics_workbench` toggle is off.
pub fn draw_acoustics_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_acoustics_workbench {
        return;
    }

    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_acoustics_workbench",
        "Acoustics",
        |app, ui| {
            ui.label(
                egui::RichText::new("native room RT60 (Sabine/Eyring) + SPL · valenx-acoustics")
                    .weak()
                    .small(),
            );
            ui.separator();

            let s = &mut app.acoustics;
            egui::ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.label(egui::RichText::new("Room").strong());
                    // Associate each numeric `DragValue` with its caption via
                    // `labelled_by`, so the spin button carries the caption as
                    // its accessibility / UI-Automation Name (egui clears a
                    // DragValue's own Name, leaving it anonymous to a screen
                    // reader / AI driver otherwise).
                    ui.horizontal(|ui| {
                        let lx = ui.label("Lx (m)");
                        ui.add(egui::DragValue::new(&mut s.length_x_m).speed(0.1))
                            .labelled_by(lx.id)
                            .on_hover_text("Room length Lx (m)");
                    });
                    ui.horizontal(|ui| {
                        let ly = ui.label("Ly (m)");
                        ui.add(egui::DragValue::new(&mut s.length_y_m).speed(0.1))
                            .labelled_by(ly.id)
                            .on_hover_text("Room width Ly (m)");
                    });
                    ui.horizontal(|ui| {
                        let lz = ui.label("Lz height (m)");
                        ui.add(egui::DragValue::new(&mut s.length_z_m).speed(0.1))
                            .labelled_by(lz.id)
                            .on_hover_text("Room height Lz (m)");
                    });
                    ui.horizontal(|ui| {
                        let ab = ui.label("absorption ᾱ");
                        ui.add(egui::DragValue::new(&mut s.absorption).speed(0.005))
                            .labelled_by(ab.id)
                            .on_hover_text("Average surface absorption coefficient ᾱ");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Air").strong());
                    ui.horizontal(|ui| {
                        let tc = ui.label("temperature (°C)");
                        ui.add(egui::DragValue::new(&mut s.temperature_c).speed(0.5))
                            .labelled_by(tc.id)
                            .on_hover_text("Air temperature (°C) — sets the speed of sound");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Source / listener").strong());
                    ui.horizontal(|ui| {
                        let sd = ui.label("source level (dB)");
                        ui.add(egui::DragValue::new(&mut s.source_db).speed(0.5))
                            .labelled_by(sd.id)
                            .on_hover_text("Source sound level (dB) at the reference distance");
                    });
                    ui.horizontal(|ui| {
                        let rd = ui.label("at distance (m)");
                        ui.add(egui::DragValue::new(&mut s.ref_distance_m).speed(0.1))
                            .labelled_by(rd.id)
                            .on_hover_text("Reference distance at which the source level is given (m)");
                    });
                    ui.horizontal(|ui| {
                        let li = ui.label("listener distance (m)");
                        ui.add(egui::DragValue::new(&mut s.listener_distance_m).speed(0.1))
                            .labelled_by(li.id)
                            .on_hover_text("Listener distance from the source (m)");
                    });

                    ui.add_space(4.0);
                    ui.label(egui::RichText::new("Monopole radiation").strong());
                    ui.label(
                        egui::RichText::new(
                            "pulsating-sphere source: peak radiated pressure + SPL at the observer",
                        )
                        .weak()
                        .small(),
                    );
                    ui.horizontal(|ui| {
                        let ra = ui.label("source radius (m)");
                        ui.add(egui::DragValue::new(&mut s.source_radius_m).speed(0.005))
                            .labelled_by(ra.id)
                            .on_hover_text("Pulsating-sphere source radius a (m)");
                    });
                    ui.horizontal(|ui| {
                        let uv = ui.label("surface velocity (m/s)");
                        ui.add(egui::DragValue::new(&mut s.surface_velocity).speed(0.001))
                            .labelled_by(uv.id)
                            .on_hover_text("Sphere surface (radial) velocity amplitude U (m/s)");
                    });
                    ui.horizontal(|ui| {
                        let fr = ui.label("frequency (Hz)");
                        ui.add(egui::DragValue::new(&mut s.frequency_hz).speed(10.0))
                            .labelled_by(fr.id)
                            .on_hover_text("Radiation frequency f (Hz)");
                    });
                    ui.horizontal(|ui| {
                        let ob = ui.label("observer distance (m)");
                        ui.add(egui::DragValue::new(&mut s.observer_distance_m).speed(0.1))
                            .labelled_by(ob.id)
                            .on_hover_text("Observer (field-point) distance r from the source (m)");
                    });

                    ui.add_space(6.0);
                    if ui
                        .button(egui::RichText::new("▶ Analyze").strong())
                        .clicked()
                    {
                        run_acoustics(s);
                    }
                    if ui
                        .button(egui::RichText::new("▶ Radiate").strong())
                        .on_hover_text(
                            "Compute the peak radiated pressure and SPL of the pulsating-sphere (monopole) source at the observer distance",
                        )
                        .clicked()
                    {
                        compute_radiation_now(s);
                    }

                    if let Some((pressure, level)) = s.radiation {
                        ui.separator();
                        ui.label(egui::RichText::new("Monopole radiation").strong());
                        ui.label(
                            egui::RichText::new(format!(
                                "radiated pressure : {pressure:.5e} Pa\nSPL                : {level:.2} dB"
                            ))
                            .monospace()
                            .small()
                            .color(egui::Color32::from_rgb(150, 200, 230)),
                        );
                    }
                    if ui
                        .button(egui::RichText::new("▶ Show 3-D room").strong())
                        .on_hover_text(
                            "Build a representative rectangular room shell with the point source inside as a 3-D solid and load it into the central viewport to orbit",
                        )
                        .clicked()
                    {
                        s.show_3d_request = true;
                    }

                    if let Some(e) = &s.error {
                        ui.add_space(4.0);
                        ui.colored_label(egui::Color32::from_rgb(220, 90, 90), e);
                    }

                    if !s.result.is_empty() {
                        ui.separator();
                        ui.label(egui::RichText::new("Room acoustics").strong());
                        ui.label(egui::RichText::new(&s.result).monospace().small());
                    }
                });
        },
    );
    if close {
        app.show_acoustics_workbench = false;
    }

    // Serviced after the panel draws (the `&mut app.acoustics` borrow is
    // released here): build the room's 3-D solid and load it.
    if app.acoustics.show_3d_request {
        app.acoustics.show_3d_request = false;
        load_room_3d(app);
    }
}

/// Validate the form, evaluate the room and format the readout.
fn run_acoustics(s: &mut AcousticsWorkbenchState) {
    match compute(s) {
        Ok(r) => {
            s.result = r;
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// Compute the **monopole radiation** answer (peak radiated pressure +
/// SPL) for the current source-sphere inputs, via
/// [`valenx_acoustics::radiation::monopole_pressure`] and
/// [`valenx_acoustics::spl`], and store it (or surface an error). Factored
/// out so the in-panel **Radiate** button and the `acoustics.compute`
/// bridge id share one path. The speed of sound follows the form's air
/// temperature; air density is the crate's [`RHO_AIR`].
fn compute_radiation_now(s: &mut AcousticsWorkbenchState) {
    match compute_radiation(s) {
        Ok(pair) => {
            s.radiation = Some(pair);
            s.error = None;
        }
        Err(e) => s.error = Some(e),
    }
}

/// The fallible monopole-radiation pipeline: peak radiated pressure `|p(r)|`
/// (Pa) and its SPL (dB re 20 µPa). Separated so it is unit-testable.
fn compute_radiation(s: &AcousticsWorkbenchState) -> Result<(f64, f64), String> {
    let c = speed_of_sound(s.temperature_c).map_err(|e| e.to_string())?;
    let pressure = monopole_pressure(
        s.source_radius_m,
        s.surface_velocity,
        s.frequency_hz,
        c,
        RHO_AIR,
        s.observer_distance_m,
    )
    .map_err(|e| e.to_string())?;
    let level = spl(pressure).map_err(|e| e.to_string())?;
    Ok((pressure, level))
}

impl AcousticsWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set
    /// via `SetControl` (returned by `ListControls`). These are the four
    /// monopole-radiation inputs the `acoustics.compute` bridge id drives —
    /// source radius + surface velocity + frequency + observer distance — so
    /// an agent can set them by caption and read back the radiated pressure
    /// + SPL.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "source radius (m)",
            "surface velocity (m/s)",
            "frequency (Hz)",
            "observer distance (m)",
        ]
    }

    /// Set one labelled control by its caption, for the agent `SetControl`
    /// bridge. Fail-loud on an unknown caption / wrong type / out-of-range;
    /// no state is written on error and nothing panics. All four radiation
    /// inputs read a positive finite `f64`.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        fn positive(v: f64, what: &str) -> Result<f64, String> {
            if !v.is_finite() || v <= 0.0 {
                return Err(format!("{what} must be finite and > 0, got {v}"));
            }
            Ok(v)
        }
        match name {
            "source radius (m)" => {
                self.source_radius_m = positive(value.as_f64()?, "source radius (m)")?
            }
            "surface velocity (m/s)" => {
                self.surface_velocity = positive(value.as_f64()?, "surface velocity (m/s)")?
            }
            "frequency (Hz)" => self.frequency_hz = positive(value.as_f64()?, "frequency (Hz)")?,
            "observer distance (m)" => {
                self.observer_distance_m = positive(value.as_f64()?, "observer distance (m)")?
            }
            other => return Err(format!("unknown acoustics control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the
    /// monopole-radiation answer (radiated pressure + SPL) once a radiation
    /// compute exists. `Some` once computed, `None` before the first
    /// `acoustics.compute` (or after a compute error with no prior result).
    pub fn agent_readout(&self) -> Option<String> {
        let (pressure, level) = self.radiation?;
        Some(format!(
            "Acoustics \u{00B7} a={:.4} m U={:.4} m/s f={:.1} Hz r={:.3} m \u{00B7} \
             p={:.5e} Pa \u{00B7} SPL={:.2} dB",
            self.source_radius_m,
            self.surface_velocity,
            self.frequency_hz,
            self.observer_distance_m,
            pressure,
            level,
        ))
    }
}

/// Run the monopole-radiation compute (the in-panel **Radiate** action).
/// Factored out so the button and the `acoustics.compute` bridge id share
/// one path.
pub(crate) fn run(app: &mut ValenxApp) {
    compute_radiation_now(&mut app.acoustics);
}

/// The room geometry `(volume V, surface area S)` in `(m³, m²)`, validated
/// through [`RoomDimensions`]. Extracted so it is unit-testable and shared
/// between the readout and the 3-D gate.
fn geometry(s: &AcousticsWorkbenchState) -> Result<(f64, f64), String> {
    RoomDimensions::new(s.length_x_m, s.length_y_m, s.length_z_m).map_err(|e| e.to_string())?;
    let v = s.length_x_m * s.length_y_m * s.length_z_m;
    let surf = 2.0
        * (s.length_x_m * s.length_y_m + s.length_y_m * s.length_z_m + s.length_z_m * s.length_x_m);
    Ok((v, surf))
}

/// Evaluate the room and format the full readout, mapping any domain error
/// to a display string. Extracted so it is unit-testable.
fn compute(s: &AcousticsWorkbenchState) -> Result<String, String> {
    let (v, surf) = geometry(s)?;
    let c = speed_of_sound(s.temperature_c).map_err(|e| e.to_string())?;

    // Total absorption A = S * a (all six surfaces share the coefficient).
    let a = total_absorption(&[(surf, s.absorption)]).map_err(|e| e.to_string())?;

    let rt_sabine = sabine_reverberation_time(v, a, c).map_err(|e| e.to_string())?;
    let rt_eyring =
        eyring_reverberation_time(v, surf, s.absorption, c).map_err(|e| e.to_string())?;

    // Lowest mode: the first entry of the ascending modal stack (the
    // axial mode along the longest wall for a plain rectangular room).
    let modes = RoomDimensions::new(s.length_x_m, s.length_y_m, s.length_z_m)
        .map_err(|e| e.to_string())?
        .modes_up_to(1, c)
        .map_err(|e| e.to_string())?;
    let low = modes.first().ok_or_else(|| "no room modes".to_string())?;
    let f_low = low.frequency_hz;
    let (nx, ny, nz) = (low.nx, low.ny, low.nz);
    let kind = match low.kind {
        valenx_acoustics::room::ModeKind::Axial => "axial",
        valenx_acoustics::room::ModeKind::Tangential => "tangential",
        valenx_acoustics::room::ModeKind::Oblique => "oblique",
    };

    // Free-field inverse-square distance drop: pressure scales as 1/r, so
    // the listener level is the source level read back after scaling its
    // RMS pressure by ref/listener distance. Reuses the verified
    // `pressure_from_spl` / `spl` pair (equivalent to the 20·log10 rule).
    if !s.ref_distance_m.is_finite() || s.ref_distance_m <= 0.0 {
        return Err("reference distance must be finite and positive".to_string());
    }
    if !s.listener_distance_m.is_finite() || s.listener_distance_m <= 0.0 {
        return Err("listener distance must be finite and positive".to_string());
    }
    let p_ref_dist = pressure_from_spl(s.source_db);
    let p_listener = p_ref_dist * (s.ref_distance_m / s.listener_distance_m);
    let level_listener = spl(p_listener).map_err(|e| e.to_string())?;
    let drop_db = s.source_db - level_listener;

    Ok(format!(
        "room Lx·Ly·Lz : {lx:.2} × {ly:.2} × {lz:.2} m\n\
         volume V       : {v:.2} m³\n\
         surface area S : {surf:.2} m²\n\
         absorption ᾱ   : {alpha:.3}\n\
         total absorp A : {a:.2} sabins\n\
         temperature    : {temp:.1} °C  (c {c:.1} m/s)\n\n\
         RT60 (Sabine)  : {rt_sabine:.2} s\n\
         RT60 (Eyring)  : {rt_eyring:.2} s\n\
         lowest mode    : {f_low:.1} Hz ({nx},{ny},{nz} {kind})\n\n\
         source level   : {src:.1} dB at {r1:.2} m\n\
         listener       : {lvl:.1} dB at {r2:.2} m\n\
         distance drop  : {drop_db:.1} dB",
        lx = s.length_x_m,
        ly = s.length_y_m,
        lz = s.length_z_m,
        alpha = s.absorption,
        temp = s.temperature_c,
        src = s.source_db,
        r1 = s.ref_distance_m,
        lvl = level_listener,
        r2 = s.listener_distance_m,
    ))
}

/// Append an outward-facing box (centre `c`, half-extents `h`) to the
/// buffers.
fn push_box(
    nodes: &mut Vec<Vector3<f64>>,
    tris: &mut Vec<usize>,
    c: Vector3<f64>,
    h: Vector3<f64>,
) {
    let base = nodes.len();
    let signs = [
        (-1.0, -1.0, -1.0),
        (1.0, -1.0, -1.0),
        (1.0, 1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, 1.0),
        (-1.0, 1.0, 1.0),
    ];
    for (sx, sy, sz) in signs {
        nodes.push(c + Vector3::new(sx * h.x, sy * h.y, sz * h.z));
    }
    let faces = [
        [1, 2, 6, 5],
        [0, 4, 7, 3],
        [3, 7, 6, 2],
        [0, 1, 5, 4],
        [4, 5, 6, 7],
        [0, 3, 2, 1],
    ];
    for f in faces {
        tris.extend_from_slice(&[
            base + f[0],
            base + f[1],
            base + f[2],
            base + f[0],
            base + f[2],
            base + f[3],
        ]);
    }
}

/// Build the room as a triangle [`Mesh`] — the rectangular room shell (a box
/// scaled to `Lx·Ly·Lz`, sitting on the floor) with a small source box
/// floating inside it. Representative geometry, true to the room's
/// proportions; the reverberation / SPL numbers are the `valenx-acoustics`
/// result. `None` for an invalid configuration.
fn room_solid_mesh(s: &AcousticsWorkbenchState) -> Option<Mesh> {
    geometry(s).ok()?;

    let mut nodes: Vec<Vector3<f64>> = Vec::new();
    let mut tris: Vec<usize> = Vec::new();

    let (lx, ly, lz) = (s.length_x_m, s.length_y_m, s.length_z_m);

    // Room shell: a box scaled to the room, resting on z = 0.
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, lz * 0.5),
        Vector3::new(lx * 0.5, ly * 0.5, lz * 0.5),
    );
    // Point source: a small box near the room centre, ~1 m off the floor.
    let src = (lx.min(ly).min(lz) * 0.08).clamp(0.05, 0.5);
    push_box(
        &mut nodes,
        &mut tris,
        Vector3::new(0.0, 0.0, lz * 0.4),
        Vector3::new(src, src, src),
    );

    let mut block = ElementBlock::new(ElementType::Tri3);
    block.connectivity = tris.iter().map(|&i| i as u32).collect();
    let mut mesh = Mesh::new("valenx-acoustics");
    mesh.nodes = nodes;
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    Some(mesh)
}

/// Build the 3-D room solid and load it into the central viewport.
fn load_room_3d(app: &mut ValenxApp) {
    let Some(mesh) = room_solid_mesh(&app.acoustics) else {
        app.acoustics.error =
            Some("room parameters are invalid — cannot build the 3-D solid".into());
        return;
    };
    let quality = valenx_mesh::quality_report(&mesh);
    let aspect_hist = valenx_mesh::aspect_ratio_histogram(&mesh, valenx_mesh::DEFAULT_AR_BUCKETS);
    let skew_hist = valenx_mesh::skewness_histogram(&mesh, valenx_mesh::DEFAULT_SKEW_BUCKETS);
    app.stl = None;
    app.mesh = Some(LoadedMesh {
        path: PathBuf::from("<room>/valenx-acoustics"),
        mesh,
        quality,
        aspect_hist,
        skew_hist,
    });
    app.frame_current_mesh();
}

/// The agent-bridge **`show_3d{kind:"acoustics"}`** product: the canonical
/// room built as a 3-D solid, paired with the workbench's own `compute()`
/// readout rows, at a fixed 3/4 camera. Registered in
/// [`crate::products_registry`]; the per-tool builder the registry dispatches
/// to. Pure — driven off [`AcousticsWorkbenchState::default`].
pub(crate) fn acoustics_product() -> crate::WorkspaceProduct {
    let s = AcousticsWorkbenchState::default();
    let mesh = room_solid_mesh(&s).expect("canonical room ⇒ solid builds");
    let loaded = crate::products_registry::loaded_mesh_from(mesh, "<room>/valenx-acoustics");
    let lines = crate::products_registry::lines_from_readout(
        &compute(&s).expect("canonical room ⇒ readout computes"),
    );
    let camera = crate::products_registry::camera_for(&loaded.mesh);
    crate::WorkspaceProduct {
        title: "Acoustics (room)".into(),
        lines,
        mesh: Some(loaded),
        vertex_colors: None,
        camera,
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
    fn default_state_is_idle() {
        let s = AcousticsWorkbenchState::default();
        assert!(s.result.is_empty());
        assert!(s.error.is_none());
    }

    #[test]
    fn analyze_default_reports_rt60_modes_and_spl() {
        let mut s = AcousticsWorkbenchState::default();
        run_acoustics(&mut s);
        assert!(
            s.error.is_none(),
            "default room should analyze: {:?}",
            s.error
        );
        assert!(s.result.contains("RT60 (Sabine)"));
        assert!(s.result.contains("RT60 (Eyring)"));
        assert!(s.result.contains("lowest mode"));
        assert!(s.result.contains("distance drop"));
        // 7x5x3 m at a=0.15, 20 C: V=105, S=142, A=21.3, Sabine ~0.79 s.
        assert!(s.result.contains("0.79"), "result: {}", s.result);
        // Eyring is shorter for this absorption: ~0.73 s.
        assert!(s.result.contains("0.73"), "result: {}", s.result);
        // 85 dB at 1 m reads ~73.0 dB at the 4 m listener.
        assert!(s.result.contains("73.0"), "result: {}", s.result);
    }

    #[test]
    fn analyze_rejects_zero_dimension() {
        let mut s = AcousticsWorkbenchState {
            length_x_m: 0.0,
            ..Default::default()
        };
        run_acoustics(&mut s);
        assert!(s.error.is_some());
    }

    #[test]
    fn sabine_matches_textbook_coefficient_ground_truth() {
        // Ground truth: the textbook RT60 = 0.161 * V / A reproduces the
        // crate's Sabine value at room temperature. For the default room
        // V = 105 m^3, A = S*a = 142*0.15 = 21.3 sabins:
        // 0.161 * 105 / 21.3 = 0.7937 s, matching `sabine_reverberation_time`.
        let s = AcousticsWorkbenchState::default();
        let (v, surf) = geometry(&s).unwrap();
        assert!((v - 105.0).abs() < 1e-9, "V = {v}");
        assert!((surf - 142.0).abs() < 1e-9, "S = {surf}");
        let a = total_absorption(&[(surf, s.absorption)]).unwrap();
        assert!((a - 21.3).abs() < 1e-9, "A = {a}");
        let c = speed_of_sound(20.0).unwrap();
        let rt = sabine_reverberation_time(v, a, c).unwrap();
        let textbook = 0.161 * v / a;
        assert!(
            (rt - textbook).abs() < 2e-3,
            "rt {rt} vs textbook {textbook}"
        );
        assert!((rt - 0.7937).abs() < 1e-3, "rt = {rt}");
    }

    #[test]
    fn spl_drop_is_six_db_per_distance_doubling_ground_truth() {
        // Ground truth: free-field inverse-square law drops 6.02 dB per
        // doubling of distance (20*log10(2)). Source at 1 m, listener at 2 m.
        let s = AcousticsWorkbenchState {
            ref_distance_m: 1.0,
            listener_distance_m: 2.0,
            ..Default::default()
        };
        let p_ref = pressure_from_spl(s.source_db);
        let p_listener = p_ref * (s.ref_distance_m / s.listener_distance_m);
        let level = spl(p_listener).unwrap();
        let drop = s.source_db - level;
        assert!((drop - 6.0206).abs() < 1e-3, "drop = {drop} dB");
    }

    /// Ground truth: the monopole (1/r) free-field law — the radiated
    /// pressure HALVES when the observer distance DOUBLES.
    #[test]
    fn monopole_pressure_halves_per_distance_doubling_ground_truth() {
        let near = AcousticsWorkbenchState {
            observer_distance_m: 1.0,
            ..Default::default()
        };
        let far = AcousticsWorkbenchState {
            observer_distance_m: 2.0,
            ..Default::default()
        };
        let (p1, _) = compute_radiation(&near).expect("near pressure");
        let (p2, _) = compute_radiation(&far).expect("far pressure");
        assert!((p2 - 0.5 * p1).abs() < 1e-9 * p1, "p1={p1} p2={p2}");
    }

    #[test]
    fn radiation_compute_produces_result_and_readout() {
        let mut s = AcousticsWorkbenchState::default();
        compute_radiation_now(&mut s);
        assert!(
            s.error.is_none(),
            "default radiation should compute: {:?}",
            s.error
        );
        let (pressure, level) = s.radiation.expect("a radiation result");
        assert!(pressure > 0.0, "pressure = {pressure}");
        assert!(level.is_finite(), "spl = {level}");
        let readout = s.agent_readout().expect("readout after radiate");
        assert!(readout.contains("Acoustics"), "readout: {readout}");
        assert!(readout.contains("SPL="), "readout names SPL: {readout}");
    }

    #[test]
    fn agent_set_radiation_controls() {
        use crate::agent_commands::AgentValue;
        let mut s = AcousticsWorkbenchState::default();
        s.agent_set("source radius (m)", &AgentValue::Float(0.2))
            .expect("a");
        assert!((s.source_radius_m - 0.2).abs() < 1e-9);
        s.agent_set("surface velocity (m/s)", &AgentValue::Float(0.05))
            .expect("U");
        assert!((s.surface_velocity - 0.05).abs() < 1e-9);
        s.agent_set("frequency (Hz)", &AgentValue::Float(500.0))
            .expect("f");
        assert!((s.frequency_hz - 500.0).abs() < 1e-9);
        s.agent_set("observer distance (m)", &AgentValue::Float(3.0))
            .expect("r");
        assert!((s.observer_distance_m - 3.0).abs() < 1e-9);
    }

    #[test]
    fn agent_set_rejects_bad_radiation_values() {
        use crate::agent_commands::AgentValue;
        let mut s = AcousticsWorkbenchState::default();
        assert!(s
            .agent_set("source radius (m)", &AgentValue::Float(-1.0))
            .is_err());
        assert!(s
            .agent_set("frequency (Hz)", &AgentValue::Float(0.0))
            .is_err());
        assert!(s
            .agent_set("observer distance (m)", &AgentValue::Float(f64::NAN))
            .is_err());
        assert!(s.agent_set("bogus", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn readout_is_none_before_radiate() {
        let s = AcousticsWorkbenchState::default();
        assert!(
            s.agent_readout().is_none(),
            "no readout before a radiation compute"
        );
    }

    #[test]
    fn run_bridge_helper_radiates_through_app() {
        let mut app = ValenxApp::default();
        run(&mut app);
        assert!(
            app.acoustics.radiation.is_some(),
            "the acoustics.compute bridge helper should produce a radiation result"
        );
    }

    #[test]
    fn control_names_are_listed() {
        let names = AcousticsWorkbenchState::agent_control_names();
        for c in [
            "source radius (m)",
            "surface velocity (m/s)",
            "frequency (Hz)",
            "observer distance (m)",
        ] {
            assert!(names.contains(&c), "missing control name {c}");
        }
    }

    #[test]
    fn room_mesh_for_default_is_nonempty_and_in_range() {
        let s = AcousticsWorkbenchState::default();
        let mesh = room_solid_mesh(&s).expect("default room yields a solid");
        assert!(mesh.nodes.len() > 8, "expected room shell + source box");
        let n = mesh.nodes.len() as u32;
        for blk in &mesh.element_blocks {
            assert!(!blk.connectivity.is_empty());
            assert_eq!(blk.connectivity.len() % 3, 0);
            assert!(blk.connectivity.iter().all(|&i| i < n));
        }
    }

    #[test]
    fn room_mesh_none_for_invalid() {
        let s = AcousticsWorkbenchState {
            length_z_m: 0.0,
            ..Default::default()
        };
        assert!(room_solid_mesh(&s).is_none());
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_workbench(app: &mut ValenxApp) {
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_acoustics_workbench(app, ctx);
        });
    }

    /// As [`draw_workbench`], but with accesskit enabled, returning the emitted
    /// accessibility tree nodes — the same tree a screen reader / AI driver
    /// consumes. `accesskit` is re-exported by egui, so no extra dependency.
    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_acoustics_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_acoustics_workbench);
        draw_workbench(&mut app);
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_acoustics_workbench = true;
        run_acoustics(&mut app.acoustics);
        draw_workbench(&mut app);
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // The room + air + source/listener DragValues are SpinButtons; each
        // must be `labelled_by` its caption (egui clears a DragValue's own
        // Name), so an AI / screen reader can find the control by the caption.
        let mut app = ValenxApp::default();
        app.show_acoustics_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // Lx, Ly, Lz, absorption, temperature, source level, ref distance,
        // listener distance.
        assert!(
            spin_buttons.len() >= 8,
            "expected the acoustics numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every acoustics DragValue must be labelled_by its caption (AI-drivable name)"
        );

        for caption in ["Lx (m)", "absorption ᾱ", "listener distance (m)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
        // The Analyze button stays a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| n.role() == Role::Button
                && n.name().is_some_and(|s| s.contains("Analyze"))),
            "the Analyze button is a named, invokable node"
        );
    }
}
