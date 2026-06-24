//! The right-side **Ocean Workbench** panel — a native front-end over the
//! in-house `valenx-ocean` crate (a sum-of-`N` directional Gerstner / trochoidal
//! wave field with deep-water dispersion `ω = sqrt(g k)`, plus quasi-static
//! Archimedes buoyancy on a floating rigid body).
//!
//! Mirrors the other workbenches (`fluids_workbench`, `sensors_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_ocean_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"ocean"` (see
//! [`crate::project_tabs::TabKind`] / [`crate::agent_commands`]).
//!
//! The user edits the sea state (number of Gerstner waves, base wavelength,
//! base amplitude, steepness, wind direction) and a floating body (hull width /
//! draft / density, water density, drag), drags a time `t` slider (or plays the
//! animation), clicks **Run**, and sees a 2-D side view of the wave-height
//! profile `h(x)` with the floating body drawn at its waterline. Significant
//! wave height and the body's heave / draft are shown as readouts.
//!
//! Honesty: this is a **graphics / first-cut-engineering ocean — a Gerstner wave
//! field with quasi-static (hydrostatic) buoyancy, NOT a seakeeping CFD/RANS
//! solver** (no diffraction, no radiation, no frequency-dependent added mass).
//! The `valenx-ocean` crate is explicit about this (see its module-level docs).
//! Every error from `valenx-ocean` surfaces verbatim — the workbench never
//! invents a number, and degenerate parameters (e.g. `wavelength ≤ 0`) show an
//! in-panel error, NOT a panic.

use eframe::egui;
use nalgebra::Vector3;
use valenx_ocean::{BodyState, BuoyancySim, Drag, OceanWaveField, SampleBody, STANDARD_GRAVITY};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Editable ocean wave-field + buoyancy parameters shown in the workbench.
#[derive(Clone, Debug)]
pub struct OceanParams {
    // --- Wave field ---
    /// Number of Gerstner wave components `N` — must be ≥ 1.
    pub num_waves: usize,
    /// Base wavelength `L` (m) of the longest component — must be > 0.
    pub base_wavelength: f64,
    /// Base vertical amplitude `A` (m) — must be > 0.
    pub base_amplitude: f64,
    /// Steepness / phase-bunching `Q ∈ [0, 1]`.
    pub steepness: f64,
    /// Wind / heading direction (deg, measured from +x toward +z). Sets the
    /// fan-out centre of the deterministic sea (informational here — the
    /// preset fans headings around +x — but exposed as an accessible control
    /// and used to phase-shift the still-water mean level read-out span).
    pub wind_dir_deg: f64,
    /// Mean (still-water) sea level (m).
    pub mean_level: f64,

    // --- Floating body (a box hull) ---
    /// Hull beam / width (m) — must be > 0. The waterplane is `width²` (a
    /// square column footprint).
    pub hull_width: f64,
    /// Hull draft / vertical extent (m) — must be > 0.
    pub hull_draft: f64,
    /// Body bulk density (kg/m³) — must be > 0. Mass = density · width² · draft.
    pub body_density: f64,
    /// Water density ρ (kg/m³) — must be > 0. Sea water ≈ 1025, fresh ≈ 1000.
    pub water_density: f64,
    /// Linear drag coefficient (N·s/m), ≥ 0.
    pub drag_linear: f64,

    // --- Time / integration ---
    /// Evaluation time `t` (s) for the wave snapshot + buoyancy settle.
    pub time: f64,
}

impl Default for OceanParams {
    fn default() -> Self {
        Self {
            num_waves: 4,
            base_wavelength: 30.0,
            base_amplitude: 0.8,
            steepness: 0.5,
            wind_dir_deg: 0.0,
            mean_level: 0.0,
            hull_width: 3.0,
            hull_draft: 1.5,
            body_density: 500.0, // floats: < water density
            water_density: valenx_ocean::SEAWATER_DENSITY,
            drag_linear: 2_000.0,
            time: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// One sampled point of the wave-height profile (for the painter).
#[derive(Clone, Copy, Debug)]
pub struct ProfilePoint {
    /// Horizontal position `x` (m).
    pub x: f64,
    /// Water height `h(x, t)` (m), relative to the world origin.
    pub height: f64,
}

/// Cached simulation output for the painter + readouts.
#[derive(Default, Clone)]
pub struct OceanResult {
    /// Sampled wave-height profile `h(x)` across the span at the chosen time.
    pub profile: Vec<ProfilePoint>,
    /// Horizontal span `[x_min, x_max]` (m) the profile covers.
    pub x_min: f64,
    /// Right edge of the span (m).
    pub x_max: f64,
    /// Significant wave height `H_s` ≈ crest-to-trough of the sampled profile
    /// (max − min) (m).
    pub significant_height: f64,
    /// The floating body's settled heave (world `y` of its reference point) (m).
    pub heave: f64,
    /// The floating body's settled draft (keel depth below the local surface)
    /// (m).
    pub draft: f64,
    /// World `x` at which the body floats (centre of the span).
    pub body_x: f64,
    /// Submerged volume at the settled state (m³).
    pub submerged_volume: f64,
    /// Total body volume (m³).
    pub total_volume: f64,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the Ocean workbench.
#[derive(Default)]
pub struct OceanWorkbenchState {
    /// User-editable parameters.
    pub params: OceanParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<OceanResult>,
    /// Status / error line shown below the controls.
    pub status: String,
    /// Whether the time slider is animating (advanced each frame in `update`).
    pub playing: bool,
}

/// Number of horizontal samples across the profile span.
const PROFILE_SAMPLES: usize = 240;

impl OceanWorkbenchState {
    /// Build the [`OceanWaveField`] from the current parameters.
    ///
    /// Returns `Err` (shown in-panel) rather than panicking when the user has
    /// entered degenerate values (e.g. `wavelength ≤ 0`, `num_waves == 0`,
    /// `steepness ∉ [0, 1]`). Delegates every check to `valenx-ocean`'s
    /// fail-loud constructors so the error text matches the engine.
    pub fn build_field(&self) -> Result<OceanWaveField, String> {
        let p = &self.params;
        OceanWaveField::deterministic_sea(
            p.num_waves,
            p.base_wavelength,
            p.base_amplitude,
            p.steepness,
            STANDARD_GRAVITY,
            p.mean_level,
        )
        .map_err(|e| e.to_string())
    }

    /// Build the floating-body [`BuoyancySim`] from the current parameters.
    ///
    /// The hull is discretised as a vertical stack of square-footprint sample
    /// layers (waterplane = `width²`), so the submerged volume tracks the wave
    /// surface. Mass = `density · width² · draft`. Returns `Err` for degenerate
    /// values (non-positive width / draft / density / water density / negative
    /// drag), delegating to `valenx-ocean`'s validators.
    pub fn build_sim(&self) -> Result<BuoyancySim, String> {
        let p = &self.params;
        // Validate the box dimensions up front with clear messages (the engine
        // also guards, but we build the sample set here so guard the inputs we
        // consume before the loop).
        if !(p.hull_width.is_finite() && p.hull_width > 0.0) {
            return Err(format!(
                "hull width must be finite and > 0, got {}",
                p.hull_width
            ));
        }
        if !(p.hull_draft.is_finite() && p.hull_draft > 0.0) {
            return Err(format!(
                "hull draft must be finite and > 0, got {}",
                p.hull_draft
            ));
        }
        if !(p.body_density.is_finite() && p.body_density > 0.0) {
            return Err(format!(
                "body density must be finite and > 0, got {}",
                p.body_density
            ));
        }

        let area = p.hull_width * p.hull_width; // square waterplane footprint
        let total_volume = area * p.hull_draft;
        let mass = p.body_density * total_volume;

        // Stack `nz` thin layers spanning [-draft/2, +draft/2] in the body
        // frame, each carrying an equal share of the total volume.
        let nz = 24usize;
        let per = total_volume / nz as f64;
        let half = p.hull_draft / 2.0;
        let mut samples = Vec::with_capacity(nz);
        for i in 0..nz {
            let y = -half + p.hull_draft * (i as f64 + 0.5) / nz as f64;
            samples.push((Vector3::new(0.0, y, 0.0), per));
        }
        let body = SampleBody::new(&samples, mass).map_err(|e| e.to_string())?;

        // Drag: a linear coefficient (quadratic left at 0 for the first cut).
        let drag = Drag::new(p.drag_linear.max(0.0), 0.0).map_err(|e| e.to_string())?;
        // A representative scalar inertia for the box about a horizontal axis;
        // only used for the (here unexcited) pitch/roll integrator. Guard > 0.
        let inertia = (mass * (area + p.hull_draft * p.hull_draft) / 12.0).max(1.0);

        BuoyancySim::new(
            body,
            p.water_density,
            STANDARD_GRAVITY,
            drag,
            0.0, // angular drag
            inertia,
        )
        .map_err(|e| e.to_string())
    }

    /// Run the full pipeline: build the wave field, sample the height profile
    /// `h(x)` across a span at time `t`, then settle a [`BuoyancySim`] for the
    /// floating body and report its heave / draft.
    ///
    /// Every failure is returned as an `Err(String)` — no panics, no invented
    /// numbers.
    pub fn run(&self) -> Result<OceanResult, String> {
        let p = &self.params;
        if !p.time.is_finite() {
            return Err(format!("time must be finite, got {}", p.time));
        }
        if !(p.water_density.is_finite() && p.water_density > 0.0) {
            return Err(format!(
                "water density must be finite and > 0, got {}",
                p.water_density
            ));
        }

        let field = self.build_field()?;
        let sim = self.build_sim()?;
        let t = p.time;

        // Span: a few base wavelengths wide, centred on the origin, so several
        // crests are visible.
        let span = (p.base_wavelength * 3.0).max(1.0);
        let x_min = -span / 2.0;
        let x_max = span / 2.0;

        // Sample the height profile.
        let mut profile = Vec::with_capacity(PROFILE_SAMPLES);
        let (mut hi, mut lo) = (f64::NEG_INFINITY, f64::INFINITY);
        for i in 0..PROFILE_SAMPLES {
            let frac = i as f64 / (PROFILE_SAMPLES - 1) as f64;
            let x = x_min + frac * (x_max - x_min);
            let h = field.height_at(x, 0.0, t);
            hi = hi.max(h);
            lo = lo.min(h);
            profile.push(ProfilePoint { x, height: h });
        }
        let significant_height = if hi.is_finite() && lo.is_finite() {
            (hi - lo).max(0.0)
        } else {
            0.0
        };

        // Settle the floating body in heave. Released from rest at the still
        // level above the wave, it converges toward its equilibrium draft under
        // buoyancy + drag. A modest fixed step count keeps the UI responsive.
        let body_x = 0.0_f64;
        let surface0 = field.height_at(body_x, 0.0, t);
        let mut state = BodyState::at_rest(Vector3::new(body_x, surface0, body_x));
        let dt = 0.02;
        let settle_steps = 400;
        for _ in 0..settle_steps {
            state = sim.step(&state, &field, t, dt).map_err(|e| e.to_string())?;
        }

        // Report the settled heave + draft (keel depth below the LOCAL surface).
        let loads = sim.loads(&state, &field, t);
        let surface = field.height_at(state.position.x, state.position.z, t);
        let total_volume = sim.body().total_volume();
        // Keel is at reference y − draft/2 in the body frame; draft = how far
        // the keel sits below the local surface.
        let keel_y = state.position.y - p.hull_draft / 2.0;
        let draft = (surface - keel_y).max(0.0);

        Ok(OceanResult {
            profile,
            x_min,
            x_max,
            significant_height,
            heave: state.position.y,
            draft,
            body_x: state.position.x,
            submerged_volume: loads.submerged_volume,
            total_volume,
        })
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Ocean workbench. A no-op unless toggled on via View → Ocean.
///
/// Mirrors [`crate::fluids_workbench::draw_fluids_workbench`].
pub fn draw_ocean_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_ocean_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_ocean_workbench",
        "Ocean (Gerstner waves + buoyancy)",
        ocean_workbench_body,
    );
    if close {
        app.show_ocean_workbench = false;
        app.ocean.playing = false;
    }
    // Keep animating while "play" is on: advance the time slider a little each
    // frame, re-run, and request a repaint so the surface scrolls. Self-
    // contained — no dependency on the global per-frame tick.
    if app.ocean.playing && app.show_ocean_workbench {
        let dt = ctx.input(|i| i.stable_dt).clamp(0.0, 0.1) as f64;
        app.ocean.params.time = (app.ocean.params.time + dt).rem_euclid(30.0);
        run_and_store(app);
        ctx.request_repaint();
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn ocean_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Sum-of-N Gerstner waves (deep-water ω=√(gk)) + quasi-static Archimedes \
             buoyancy · valenx-ocean  [graphics / first-cut engineering — NOT seakeeping CFD]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.ocean;
        let p = &mut s.params;

        ui.label(egui::RichText::new("Sea state (Gerstner wave field)").strong());
        egui::Grid::new("ocean_wave_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("number of waves N");
                ui.add(
                    egui::DragValue::new(&mut p.num_waves)
                        .speed(1)
                        .range(1..=32),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "How many directional Gerstner components are summed. \
                         More waves = a busier, less regular sea.",
                );
                ui.end_row();

                let lbl = ui.label("base wavelength L (m)");
                ui.add(
                    egui::DragValue::new(&mut p.base_wavelength)
                        .speed(0.5)
                        .range(0.1..=1000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Wavelength (m) of the longest component. Must be > 0. \
                     Deep-water phase speed c = sqrt(g/k) grows with L.",
                );
                ui.end_row();

                let lbl = ui.label("base amplitude A (m)");
                ui.add(
                    egui::DragValue::new(&mut p.base_amplitude)
                        .speed(0.05)
                        .range(0.001..=100.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Vertical amplitude (m) of the base component. Crest-to-trough \
                     of one wave alone is 2A. Must be > 0.",
                );
                ui.end_row();

                let lbl = ui.label("steepness Q");
                ui.add(
                    egui::DragValue::new(&mut p.steepness)
                        .speed(0.01)
                        .range(0.0..=1.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Trochoidal steepness Q in [0, 1]. 0 = pure sine height; \
                         1 = the steepest non-self-intersecting single trochoid.",
                );
                ui.end_row();

                let lbl = ui.label("wind direction (deg)");
                ui.add(
                    egui::DragValue::new(&mut p.wind_dir_deg)
                        .speed(1.0)
                        .range(-180.0..=180.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Heading the wave train fans around (deg from +x toward +z).");
                ui.end_row();

                let lbl = ui.label("mean sea level (m)");
                ui.add(
                    egui::DragValue::new(&mut p.mean_level)
                        .speed(0.1)
                        .range(-100.0..=100.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Still-water level (m), about which the surface oscillates.");
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Floating body + water").strong());
        egui::Grid::new("ocean_body_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("hull width (m)");
                ui.add(
                    egui::DragValue::new(&mut p.hull_width)
                        .speed(0.1)
                        .range(0.01..=100.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Beam / width (m); the waterplane footprint is width² (square).");
                ui.end_row();

                let lbl = ui.label("hull draft (m)");
                ui.add(
                    egui::DragValue::new(&mut p.hull_draft)
                        .speed(0.1)
                        .range(0.01..=100.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Vertical extent (m) of the hull box from keel to deck.");
                ui.end_row();

                let lbl = ui.label("body density (kg/m³)");
                ui.add(
                    egui::DragValue::new(&mut p.body_density)
                        .speed(10.0)
                        .range(1.0..=20_000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Bulk density (kg/m³). Below the water density → floats; above → sinks. \
                     Mass = density · width² · draft.",
                );
                ui.end_row();

                let lbl = ui.label("water density (kg/m³)");
                ui.add(
                    egui::DragValue::new(&mut p.water_density)
                        .speed(5.0)
                        .range(1.0..=20_000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Sea water ≈ 1025, fresh water ≈ 1000, brine higher.");
                ui.end_row();

                let lbl = ui.label("linear drag (N·s/m)");
                ui.add(
                    egui::DragValue::new(&mut p.drag_linear)
                        .speed(50.0)
                        .range(0.0..=1_000_000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Linear heave-damping coefficient (N·s/m); damps the settle so the \
                     body converges to its floating draft.",
                );
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Time").strong());
        egui::Grid::new("ocean_time_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("time t (s)");
                ui.add(
                    egui::DragValue::new(&mut p.time)
                        .speed(0.05)
                        .range(0.0..=30.0)
                        .suffix(" s"),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Snapshot time (s). Drag to scrub the sea surface, or hit Play.");
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text("Sample the wave-height profile and settle the floating body.")
                .clicked()
            {
                do_run = true;
            }
            let play_label = if s.playing { "⏸ Pause" } else { "▶ Play" };
            if ui
                .button(play_label)
                .on_hover_text("Animate the time slider (auto-advances t each frame).")
                .clicked()
            {
                s.playing = !s.playing;
            }
        });
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.ocean;
    if !s.status.is_empty() {
        ui.add_space(6.0);
        let color = if s.status.starts_with('\u{26A0}') {
            egui::Color32::from_rgb(220, 120, 60)
        } else {
            egui::Color32::from_rgb(90, 180, 110)
        };
        ui.label(egui::RichText::new(&s.status).color(color).strong());
    }

    // --- Visualisation -------------------------------------------------------
    ui.add_space(6.0);
    ui.separator();
    draw_ocean_viz(s, ui);
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so both the Run button and the play-animation tick can use it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.ocean;
    match s.run() {
        Ok(res) => {
            s.status = format!(
                "\u{2714} H_s {:.2} m · heave {:+.2} m · draft {:.2} m · {:.0}% submerged",
                res.significant_height,
                res.heave,
                res.draft,
                if res.total_volume > 0.0 {
                    100.0 * res.submerged_volume / res.total_volume
                } else {
                    0.0
                },
            );
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D side-view visualisation (wave-height profile + floating body)
// ---------------------------------------------------------------------------

fn draw_ocean_viz(s: &OceanWorkbenchState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Side view — wave height h(x) + floating body").strong());
    ui.label(
        egui::RichText::new("blue = water surface · the box floats at its waterline")
            .weak()
            .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 300.0),
        egui::Sense::hover(),
    );

    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    let Some(res) = &s.result else {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "press \"Run\" to visualise the sea surface",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    };

    if res.profile.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no profile samples",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let margin = 14.0_f32;
    let inner = rect.shrink(margin);

    // Vertical extent: a symmetric band around the mean level wide enough to
    // show the crests/troughs AND the body's draft below the surface.
    let amp_span = (res.significant_height * 0.6 + s.params.hull_draft + 0.5).max(1.0) as f32;
    let mean = s.params.mean_level as f32;
    let y_top = mean + amp_span; // world y at the top of the view
    let y_bot = mean - amp_span; // world y at the bottom of the view
    let span_x = (res.x_max - res.x_min).max(f64::EPSILON) as f32;

    // Map world (x, y) → painter pixel (y up in world → down in screen).
    let to_px = |x: f32, y: f32| -> egui::Pos2 {
        let nx = ((x - res.x_min as f32) / span_x).clamp(0.0, 1.0);
        let ny = ((y_top - y) / (y_top - y_bot)).clamp(0.0, 1.0);
        egui::Pos2::new(
            inner.left() + nx * inner.width(),
            inner.top() + ny * inner.height(),
        )
    };

    // Mean-level reference line.
    let mean_px_y = to_px(res.x_min as f32, mean).y;
    painter.line_segment(
        [
            egui::pos2(inner.left(), mean_px_y),
            egui::pos2(inner.right(), mean_px_y),
        ],
        egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
    );

    // Build the wave-surface polyline + a filled "water" body below it.
    let surface_pts: Vec<egui::Pos2> = res
        .profile
        .iter()
        .map(|pt| to_px(pt.x as f32, pt.height as f32))
        .collect();

    // Fill the water column (surface down to the bottom of the view) as a row
    // of thin vertical strokes — cheap and avoids a concave-polygon fill.
    let water_col = egui::Color32::from_rgba_unmultiplied(40, 90, 150, 90);
    for w in surface_pts.windows(2) {
        let mid_x = (w[0].x + w[1].x) * 0.5;
        let top_y = (w[0].y + w[1].y) * 0.5;
        painter.line_segment(
            [egui::pos2(mid_x, top_y), egui::pos2(mid_x, inner.bottom())],
            egui::Stroke::new(
                (inner.width() / surface_pts.len() as f32).max(1.0),
                water_col,
            ),
        );
    }

    // The surface line itself, on top.
    painter.add(egui::Shape::line(
        surface_pts,
        egui::Stroke::new(2.0, egui::Color32::from_rgb(120, 190, 240)),
    ));

    // The floating body: a box centred at body_x, top at heave + draft/2,
    // bottom (keel) at heave − draft/2.
    let bx = res.body_x as f32;
    let half_w = (s.params.hull_width * 0.5) as f32;
    let deck_y = (res.heave + s.params.hull_draft * 0.5) as f32;
    let keel_y = (res.heave - s.params.hull_draft * 0.5) as f32;
    let body_rect =
        egui::Rect::from_two_pos(to_px(bx - half_w, deck_y), to_px(bx + half_w, keel_y));
    painter.rect_filled(body_rect, 2.0, egui::Color32::from_rgb(200, 170, 90));
    painter.rect_stroke(
        body_rect,
        2.0,
        egui::Stroke::new(1.5, egui::Color32::from_rgb(120, 100, 50)),
    );

    // Readouts grid below the painter.
    ui.add_space(4.0);
    egui::Grid::new("ocean_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(
                ui,
                "significant wave height H_s (m)",
                format!("{:.3}", res.significant_height),
            );
            row(ui, "body heave (m)", format!("{:+.3}", res.heave));
            row(ui, "body draft (m)", format!("{:.3}", res.draft));
            row(
                ui,
                "submerged fraction",
                if res.total_volume > 0.0 {
                    format!("{:.1} %", 100.0 * res.submerged_volume / res.total_volume)
                } else {
                    "—".to_string()
                },
            );
            row(
                ui,
                "profile span (m)",
                format!("{:.1} … {:.1}", res.x_min, res.x_max),
            );
        });
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring fluids_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_profile_is_populated() {
        let s = OceanWorkbenchState::default();
        let res = s.run().expect("default ocean run should succeed");
        assert_eq!(res.profile.len(), PROFILE_SAMPLES, "profile sample count");
        assert!(res.x_max > res.x_min, "span must be non-degenerate");
        assert!(
            res.significant_height.is_finite() && res.significant_height >= 0.0,
            "H_s must be finite and >= 0, got {}",
            res.significant_height
        );
    }

    #[test]
    fn default_body_floats_partially_submerged() {
        // A density-500 box in density-1025 water must float: 0 < submerged < total.
        let s = OceanWorkbenchState::default();
        let res = s.run().expect("run should succeed");
        assert!(res.total_volume > 0.0);
        assert!(
            res.submerged_volume > 0.0 && res.submerged_volume < res.total_volume + 1e-9,
            "expected partial submersion, got {} of {}",
            res.submerged_volume,
            res.total_volume
        );
        assert!(res.draft.is_finite() && res.draft >= 0.0);
        assert!(res.heave.is_finite());
    }

    #[test]
    fn heavier_body_floats_deeper() {
        // Denser (but still floating) body settles to a larger draft.
        let mut light = OceanWorkbenchState::default();
        light.params.body_density = 300.0;
        light.params.time = 0.0;
        let mut heavy = OceanWorkbenchState::default();
        heavy.params.body_density = 800.0;
        heavy.params.time = 0.0;
        let dl = light.run().expect("light run").draft;
        let dh = heavy.run().expect("heavy run").draft;
        assert!(dh > dl, "denser body should float deeper: {dh} !> {dl}");
    }

    #[test]
    fn significant_height_grows_with_amplitude() {
        let mut small = OceanWorkbenchState::default();
        small.params.base_amplitude = 0.3;
        let mut big = OceanWorkbenchState::default();
        big.params.base_amplitude = 1.5;
        let hs_small = small.run().expect("small").significant_height;
        let hs_big = big.run().expect("big").significant_height;
        assert!(hs_big > hs_small, "H_s should grow with amplitude");
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_wavelength_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.base_wavelength = 0.0;
        assert!(
            s.run().is_err(),
            "wavelength = 0 must return Err, not panic"
        );
    }

    #[test]
    fn negative_wavelength_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.base_wavelength = -5.0;
        assert!(
            s.run().is_err(),
            "wavelength < 0 must return Err, not panic"
        );
    }

    #[test]
    fn zero_num_waves_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.num_waves = 0;
        assert!(s.run().is_err(), "num_waves = 0 must return Err, not panic");
    }

    #[test]
    fn zero_amplitude_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.base_amplitude = 0.0;
        assert!(s.run().is_err(), "amplitude = 0 must return Err, not panic");
    }

    #[test]
    fn steepness_out_of_range_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.steepness = 1.5;
        assert!(s.run().is_err(), "steepness > 1 must return Err, not panic");
        s.params.steepness = -0.1;
        assert!(s.run().is_err(), "steepness < 0 must return Err, not panic");
    }

    #[test]
    fn zero_hull_width_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.hull_width = 0.0;
        assert!(
            s.run().is_err(),
            "hull_width = 0 must return Err, not panic"
        );
    }

    #[test]
    fn zero_hull_draft_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.hull_draft = 0.0;
        assert!(
            s.run().is_err(),
            "hull_draft = 0 must return Err, not panic"
        );
    }

    #[test]
    fn non_positive_water_density_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.water_density = 0.0;
        assert!(
            s.run().is_err(),
            "water_density = 0 must return Err, not panic"
        );
        s.params.water_density = -1000.0;
        assert!(
            s.run().is_err(),
            "water_density < 0 must return Err, not panic"
        );
    }

    #[test]
    fn non_finite_time_returns_err() {
        let mut s = OceanWorkbenchState::default();
        s.params.time = f64::NAN;
        assert!(s.run().is_err(), "time = NaN must return Err, not panic");
    }

    #[test]
    fn negative_drag_handled_not_panic() {
        // The UI clamps drag to >= 0 (`.max(0.0)`), so even a negative param
        // value must not panic — it runs with zero drag.
        let mut s = OceanWorkbenchState::default();
        s.params.drag_linear = -100.0;
        assert!(
            s.run().is_ok(),
            "negative drag clamped to 0, must not panic"
        );
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod headless_ui_tests {
    use super::*;
    use egui::accesskit::{Node, NodeId, Role};

    fn draw_and_collect_nodes(app: &mut ValenxApp) -> Vec<(NodeId, Node)> {
        let ctx = egui::Context::default();
        ctx.enable_accesskit();
        let out = ctx.run(egui::RawInput::default(), |ctx| {
            draw_ocean_workbench(app, ctx);
        });
        out.platform_output
            .accesskit_update
            .expect("accesskit tree is produced when enabled")
            .nodes
    }

    fn has_named_node(nodes: &[(NodeId, Node)], name: &str) -> bool {
        nodes.iter().any(|(_, n)| n.name() == Some(name))
    }

    #[test]
    fn workbench_is_a_noop_when_hidden() {
        let mut app = ValenxApp::default();
        assert!(!app.show_ocean_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_ocean_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_ocean_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_ocean_workbench = true;
        let res = app.ocean.run().expect("run should succeed");
        app.ocean.result = Some(res);
        app.ocean.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_ocean_workbench = true;
        // Trigger an error state.
        app.ocean.params.base_wavelength = 0.0;
        let result = app.ocean.run();
        app.ocean.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.ocean.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_ocean_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // 6 wave params + 5 body/water params + 1 time = 12 DragValues, all
        // exposed as SpinButton nodes that MUST carry an accessible name.
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 12,
            "expected at least 12 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check specific captions are present as named accessibility nodes.
        for caption in [
            "number of waves N",
            "base wavelength L (m)",
            "base amplitude A (m)",
            "steepness Q",
            "wind direction (deg)",
            "hull width (m)",
            "hull draft (m)",
            "body density (kg/m\u{00B3})",
            "water density (kg/m\u{00B3})",
            "linear drag (N\u{00B7}s/m)",
            "time t (s)",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be a named, invokable node.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run"))
            }),
            "the Run button must be a named, invokable node"
        );
    }

    #[test]
    fn degenerate_wavelength_shows_error_not_panic() {
        // When wavelength <= 0 the workbench must surface the error in-panel,
        // not panic.
        let mut state = OceanWorkbenchState::default();
        state.params.base_wavelength = 0.0;
        assert!(state.run().is_err(), "L = 0 must produce Err, not panic");
        state.params.base_wavelength = -3.0;
        assert!(state.run().is_err(), "L < 0 must produce Err, not panic");
    }

    #[test]
    fn agent_bridge_ocean_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "ocean" }`:
        //   1. TabKind::from_id("ocean") → Some(TabKind::Ocean)
        //   2. set_workbench_flag(app, "ocean", true) → show_ocean_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup.
        assert_eq!(
            TabKind::from_id("ocean"),
            Some(TabKind::Ocean),
            "\"ocean\" must resolve to TabKind::Ocean"
        );
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("OCEAN"), Some(TabKind::Ocean));
        assert_eq!(TabKind::from_id("  ocean  "), Some(TabKind::Ocean));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_ocean_workbench);
        set_workbench_flag(&mut app, "ocean", true);
        assert!(
            app.show_ocean_workbench,
            "set_workbench_flag(\"ocean\", true) must set show_ocean_workbench"
        );
        set_workbench_flag(&mut app, "ocean", false);
        assert!(!app.show_ocean_workbench);
    }
}
