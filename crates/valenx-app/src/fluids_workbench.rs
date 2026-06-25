//! The right-side **Fluids (SPH) Workbench** panel — a native front-end over the
//! in-house `valenx-fluids` crate (particle-based SPH fluid simulation).
//!
//! Mirrors the other workbenches (`sensors_workbench`, `rotor_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_fluids_workbench`], toggled from the View menu
//! and openable by the agent bridge under the workbench id `"fluids"` (see
//! [`crate::project_tabs::TabKind`] / [`crate::agent_commands`]).
//!
//! The user edits SPH configuration (smoothing length `h`, particle mass,
//! viscosity, gravity, rest density, timestep `dt`, step count, and initial
//! fluid-block size N×N×N), clicks **Run simulation**, and sees the particle
//! cloud visualised in a 2-D side view coloured by speed. Particle count,
//! mean density, and elapsed step count are shown as readouts.
//!
//! Honesty: the simulation is **graphics / real-time-grade interactive SPH —
//! NOT validated against analytic CFD or experiment**. The `valenx-fluids`
//! crate is explicit about this (see its module-level docs). Every error from
//! `valenx-fluids` surfaces verbatim — the workbench never invents a number,
//! and degenerate parameters (e.g. `h ≤ 0`) show an in-panel error, NOT a panic.

use eframe::egui;
use nalgebra::Vector3;
use valenx_fluids::{BoxBoundary, Particle, ParticleSystem, SphConfig, SphSolver};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Editable SPH simulation parameters shown in the workbench.
#[derive(Clone, Debug)]
pub struct FluidsParams {
    /// Smoothing length `h` (m) — must be > 0.
    pub smoothing_length: f64,
    /// Per-particle mass (kg) — must be > 0. When set to `0.0` the workbench
    /// uses the water-preset mass derived from `h`.
    pub mass_override: f64,
    /// Dynamic viscosity μ (Pa·s), ≥ 0.
    pub viscosity: f64,
    /// Gravity magnitude (m/s²); direction is −Z (down in the 2-D view).
    pub gravity_z: f64,
    /// Rest density ρ₀ (kg/m³) — must be > 0.
    pub rest_density: f64,
    /// Integration time step dt (s) — must be > 0.
    pub dt: f64,
    /// Number of solver steps to run per "Run simulation" click.
    pub num_steps: usize,
    /// Particles per axis of the initial lattice block (total = n³).
    pub n_per_axis: usize,
    /// Box half-size (m) — the boundary box spans [0, box_size]³.
    pub box_size: f64,
}

impl Default for FluidsParams {
    fn default() -> Self {
        Self {
            smoothing_length: 0.1,
            mass_override: 0.0, // 0 = use SphConfig::water preset
            viscosity: 0.1,
            gravity_z: 9.806_65,
            rest_density: 1000.0,
            dt: 1e-3,
            num_steps: 50,
            n_per_axis: 4,
            box_size: 0.5,
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// A single particle snapshot for visualisation (position + speed).
#[derive(Clone, Debug)]
pub struct ParticleSnapshot {
    /// Position (m).
    pub position: Vector3<f64>,
    /// Speed |v| (m/s).
    pub speed: f64,
    /// Smoothed density (kg/m³).
    pub density: f64,
}

/// Cached simulation output for the painter.
#[derive(Default, Clone)]
pub struct FluidsResult {
    /// Per-particle snapshots.
    pub particles: Vec<ParticleSnapshot>,
    /// Total number of steps actually executed.
    pub steps_done: usize,
    /// Mean density across all particles (kg/m³).
    pub mean_density: f64,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the Fluids (SPH) workbench.
#[derive(Default)]
pub struct FluidsWorkbenchState {
    /// User-editable simulation parameters.
    pub params: FluidsParams,
    /// Last successful simulation result (populated after a successful run).
    pub result: Option<FluidsResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl FluidsWorkbenchState {
    /// Build the [`SphConfig`] from the current parameters.
    ///
    /// Returns `Err` (shown in-panel) rather than panicking when the user
    /// has entered degenerate values (e.g. `h ≤ 0`, zero mass, negative dt).
    pub fn build_config(&self) -> Result<SphConfig, String> {
        let p = &self.params;

        // Start from the water preset for the smoothing length; this
        // validates h up front and gives a sensible mass if mass_override ≤ 0.
        let mut cfg = SphConfig::water(p.smoothing_length).map_err(|e| e.to_string())?;

        // Override individual fields if the user has set non-default values.
        if p.mass_override > 0.0 {
            if !p.mass_override.is_finite() {
                return Err(format!(
                    "particle mass must be finite and > 0, got {}",
                    p.mass_override
                ));
            }
            cfg.mass = p.mass_override;
        }
        if !p.viscosity.is_finite() || p.viscosity < 0.0 {
            return Err(format!(
                "viscosity must be ≥ 0 and finite, got {}",
                p.viscosity
            ));
        }
        cfg.viscosity = p.viscosity;

        if !p.gravity_z.is_finite() {
            return Err(format!("gravity must be finite, got {}", p.gravity_z));
        }
        cfg.gravity = Vector3::new(0.0, 0.0, -(p.gravity_z));

        // Override the EOS rest density if non-default.
        if !(p.rest_density.is_finite() && p.rest_density > 0.0) {
            return Err(format!(
                "rest density must be finite and > 0, got {}",
                p.rest_density
            ));
        }
        // Rebuild the EOS with the user-supplied rest density.
        cfg.eos = valenx_fluids::EquationOfState::tait_water(p.rest_density, 50.0)
            .map_err(|e| e.to_string())?;

        Ok(cfg)
    }

    /// Build the initial [`ParticleSystem`]: an N×N×N lattice block near the
    /// top of the box, dropping under gravity.
    ///
    /// Returns `Err` for degenerate parameters (n = 0, non-finite box_size).
    pub fn build_system(&self) -> Result<ParticleSystem, String> {
        let p = &self.params;
        if p.n_per_axis == 0 {
            return Err("n particles per axis must be ≥ 1".to_string());
        }
        if !(p.box_size.is_finite() && p.box_size > 0.0) {
            return Err(format!(
                "box size must be finite and > 0, got {}",
                p.box_size
            ));
        }

        let spacing = p.smoothing_length * 0.5;
        let n = p.n_per_axis;
        let mut sys = ParticleSystem::new();
        // Place the block in the upper part of the box with a small margin.
        let origin = Vector3::new(
            0.05 * p.box_size,
            0.05 * p.box_size,
            p.box_size * 0.5, // midpoint height — gravity pulls down (−Z)
        );
        for ix in 0..n {
            for iy in 0..n {
                for iz in 0..n {
                    let pos = origin
                        + Vector3::new(
                            ix as f64 * spacing,
                            iy as f64 * spacing,
                            iz as f64 * spacing,
                        );
                    sys.push(Particle::at(pos)).map_err(|e| e.to_string())?;
                }
            }
        }
        if sys.is_empty() {
            return Err("no particles were created — check n per axis".to_string());
        }
        Ok(sys)
    }

    /// Build the [`BoxBoundary`] from the current parameters.
    pub fn build_boundary(&self) -> Result<BoxBoundary, String> {
        let s = self.params.box_size;
        if !(s.is_finite() && s > 0.0) {
            return Err(format!("box size must be finite and > 0, got {s}"));
        }
        BoxBoundary::new(Vector3::zeros(), Vector3::new(s, s, s), 0.3).map_err(|e| e.to_string())
    }

    /// Run the full simulation pipeline (build config → system → boundary →
    /// step N times) and return a [`FluidsResult`].
    ///
    /// Every failure is returned as an `Err(String)` — no panics.
    pub fn run(&self) -> Result<FluidsResult, String> {
        let p = &self.params;
        if !(p.dt.is_finite() && p.dt > 0.0) {
            return Err(format!("dt must be finite and > 0, got {}", p.dt));
        }
        if p.num_steps == 0 {
            return Err("num_steps must be ≥ 1".to_string());
        }

        let cfg = self.build_config()?;
        let mut solver = SphSolver::new(cfg).map_err(|e| e.to_string())?;
        let mut system = self.build_system()?;
        let boundary = self.build_boundary()?;

        for _ in 0..p.num_steps {
            solver
                .step(&mut system, p.dt, Some(&boundary))
                .map_err(|e| e.to_string())?;
        }

        let particles: Vec<ParticleSnapshot> = system
            .particles()
            .iter()
            .map(|par| ParticleSnapshot {
                position: par.position,
                speed: par.velocity.norm(),
                density: par.density,
            })
            .collect();

        let mean_density = if particles.is_empty() {
            0.0
        } else {
            particles.iter().map(|s| s.density).sum::<f64>() / particles.len() as f64
        };

        Ok(FluidsResult {
            particles,
            steps_done: p.num_steps,
            mean_density,
        })
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`;
    /// each string matches exactly the caption the form draws.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "smoothing length h (m)",
            "viscosity μ (Pa·s)",
            "gravity |g| (m/s²)",
            "rest density ρ₀ (kg/m³)",
            "particle mass (kg)  [0 = auto]",
            "time step dt (s)",
            "number of steps",
            "particles per axis N",
            "box size (m)",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the
    /// wrong type / out of range returns `Err(String)` — never a panic. Ranges
    /// mirror the form's `DragValue` clamps exactly.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let ranged = |v: f64, lo: f64, hi: f64, what: &str| -> Result<f64, String> {
            if v.is_finite() && (lo..=hi).contains(&v) {
                Ok(v)
            } else {
                Err(format!("{what} must be in {lo}..={hi}, got {v}"))
            }
        };
        let ranged_int = |value: &crate::agent_commands::AgentValue,
                          lo: i64,
                          hi: i64,
                          what: &str|
         -> Result<usize, String> {
            let n = value.as_i64()?;
            if (lo..=hi).contains(&n) {
                Ok(n as usize)
            } else {
                Err(format!("{what} must be in {lo}..={hi}, got {n}"))
            }
        };
        let p = &mut self.params;
        match name {
            "smoothing length h (m)" => {
                p.smoothing_length = ranged(value.as_f64()?, 1e-4, 1.0, "smoothing length h")?
            }
            "viscosity μ (Pa·s)" => {
                p.viscosity = ranged(value.as_f64()?, 0.0, 10.0, "viscosity μ")?
            }
            "gravity |g| (m/s²)" => p.gravity_z = ranged(value.as_f64()?, 0.0, 100.0, "gravity")?,
            "rest density ρ₀ (kg/m³)" => {
                p.rest_density = ranged(value.as_f64()?, 1.0, 20_000.0, "rest density ρ₀")?
            }
            "particle mass (kg)  [0 = auto]" => {
                p.mass_override = ranged(value.as_f64()?, 0.0, 100.0, "particle mass")?
            }
            "time step dt (s)" => p.dt = ranged(value.as_f64()?, 1e-6, 0.1, "time step dt")?,
            "number of steps" => p.num_steps = ranged_int(value, 1, 500, "number of steps")?,
            "particles per axis N" => {
                p.n_per_axis = ranged_int(value, 1, 10, "particles per axis N")?
            }
            "box size (m)" => p.box_size = ranged(value.as_f64()?, 0.01, 10.0, "box size")?,
            other => return Err(format!("unknown Fluids control: {other:?}")),
        }
        Ok(())
    }

    /// The current computed-result text for the agent `ReadReadout` bridge (see
    /// [`crate::agent_commands`]). This workbench keeps its result as a structured
    /// [`FluidsResult`] and renders a one-line `status` summary (a `✔ …` line on
    /// success, a `⚠ …` line on error) — that same `status` string is returned
    /// here. `None` when it is empty, i.e. the pipeline has not been run yet.
    /// Read-only — lets an agent read the answer back after driving a run,
    /// closing the live-driving loop.
    pub fn agent_readout(&self) -> Option<String> {
        if self.status.is_empty() {
            None
        } else {
            Some(self.status.clone())
        }
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Fluids (SPH) workbench. A no-op unless toggled on via View → SPH Fluids.
///
/// Mirrors [`crate::sensors_workbench::draw_sensors_workbench`].
pub fn draw_fluids_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_fluids_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_fluids_workbench",
        "Fluids (SPH particle sim)",
        fluids_workbench_body,
    );
    if close {
        app.show_fluids_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn fluids_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "SPH particle fluid (Müller et al. 2003) · valenx-fluids  \
             [graphics-grade — NOT validated CFD]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.fluids;
        let p = &mut s.params;

        ui.label(egui::RichText::new("SPH configuration").strong());
        egui::Grid::new("fluids_sph_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("smoothing length h (m)");
                ui.add(
                    egui::DragValue::new(&mut p.smoothing_length)
                        .speed(0.005)
                        .range(1e-4..=1.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Kernel support radius h (m). Must be > 0. \
                     Smaller h = sharper features but more particles needed for coverage.",
                );
                ui.end_row();

                let lbl = ui.label("viscosity μ (Pa·s)");
                ui.add(
                    egui::DragValue::new(&mut p.viscosity)
                        .speed(0.005)
                        .range(0.0..=10.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Dynamic viscosity (Pa·s); 0 = inviscid, ~0.001 = water, 1.0 = oil.",
                );
                ui.end_row();

                let lbl = ui.label("gravity |g| (m/s²)");
                ui.add(
                    egui::DragValue::new(&mut p.gravity_z)
                        .speed(0.1)
                        .range(0.0..=100.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Gravity magnitude (m/s²); direction is −Z (downward in the 2-D view).",
                );
                ui.end_row();

                let lbl = ui.label("rest density ρ₀ (kg/m³)");
                ui.add(
                    egui::DragValue::new(&mut p.rest_density)
                        .speed(10.0)
                        .range(1.0..=20_000.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "EOS rest density (kg/m³). Water ≈ 1000, sea water ≈ 1025, mercury ≈ 13 534.",
                );
                ui.end_row();

                let lbl = ui.label("particle mass (kg)  [0 = auto]");
                ui.add(
                    egui::DragValue::new(&mut p.mass_override)
                        .speed(0.0001)
                        .range(0.0..=100.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Per-particle mass (kg). Set to 0 to use the water-preset mass \
                     derived from h (≈ ρ₀·(h/2)³).",
                );
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Integration").strong());
        egui::Grid::new("fluids_integration_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("time step dt (s)");
                ui.add(
                    egui::DragValue::new(&mut p.dt)
                        .speed(0.0001)
                        .range(1e-6..=0.1),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Integration time step (s). CFL-stable limit depends on h and sound speed. \
                     1e-3 is typical for h = 0.1 m.",
                );
                ui.end_row();

                let lbl = ui.label("number of steps");
                ui.add(
                    egui::DragValue::new(&mut p.num_steps)
                        .speed(1)
                        .range(1..=500),
                )
                .labelled_by(lbl.id)
                .on_hover_text("How many solver steps to run per click.");
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Scene").strong());
        egui::Grid::new("fluids_scene_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("particles per axis N");
                ui.add(
                    egui::DragValue::new(&mut p.n_per_axis)
                        .speed(1)
                        .range(1..=10),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Particles per axis of the initial lattice block (total = N³). \
                     N = 4 → 64 particles.",
                );
                ui.end_row();

                let lbl = ui.label("box size (m)");
                ui.add(
                    egui::DragValue::new(&mut p.box_size)
                        .speed(0.01)
                        .range(0.01..=10.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Side length of the cubic boundary box (m).");
                ui.end_row();
            });

        ui.add_space(6.0);
        let n3 = (p.n_per_axis as u64).pow(3);
        if ui
            .button(
                egui::RichText::new(format!(
                    "Run simulation  ({n3} particles × {} steps)",
                    p.num_steps
                ))
                .strong(),
            )
            .on_hover_text("Build the SPH config + particle system and step the simulation.")
            .clicked()
        {
            do_run = true;
        }
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        let s = &mut app.fluids;
        match s.run() {
            Ok(res) => {
                s.status = format!(
                    "\u{2714} {} particles · mean density {:.1} kg/m\u{00B3} · {} steps",
                    res.particles.len(),
                    res.mean_density,
                    res.steps_done,
                );
                s.result = Some(res);
            }
            Err(e) => {
                s.status = format!("\u{26A0} {e}");
                s.result = None;
            }
        }
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.fluids;
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
    draw_particle_viz(s, ui);
}

// ---------------------------------------------------------------------------
// 2-D particle visualisation (X–Z side view, coloured by speed)
// ---------------------------------------------------------------------------

fn draw_particle_viz(s: &FluidsWorkbenchState, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Particle cloud — X/Z side view (colour = speed)").strong());
    ui.label(
        egui::RichText::new("blue = slow · red = fast · dot size proportional to density")
            .weak()
            .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(420.0), 300.0),
        egui::Sense::hover(),
    );

    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(16, 18, 24));

    let Some(res) = &s.result else {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "press \"Run simulation\" to visualise",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(110),
        );
        return;
    };

    if res.particles.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no particles",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(110),
        );
        return;
    }

    // Determine scene bounds from particle positions.
    let box_s = s.params.box_size as f32;
    let margin = 12.0_f32;
    let inner = rect.shrink(margin);

    // Map (X, Z) scene coords → painter pixel.
    // X → horizontal, Z → vertical (Z = 0 at top, Z = box_size at bottom —
    // inverted so the fluid falls *downward* in the view as gravity pulls it).
    let to_px = |x: f32, z: f32| -> egui::Pos2 {
        let nx = (x / box_s.max(f32::EPSILON)).clamp(0.0, 1.0);
        let nz = 1.0 - (z / box_s.max(f32::EPSILON)).clamp(0.0, 1.0);
        egui::Pos2::new(
            inner.left() + nx * inner.width(),
            inner.top() + nz * inner.height(),
        )
    };

    // Compute max speed for colour normalisation (clamp ≥ 1e-6 to avoid div-0).
    let max_speed = res
        .particles
        .iter()
        .map(|p| p.speed)
        .fold(0.0_f64, f64::max)
        .max(1e-6) as f32;

    // Compute min/max density for dot-size scaling.
    let (min_d, max_d) = res
        .particles
        .iter()
        .fold((f64::MAX, 0.0_f64), |(mn, mx), p| {
            (mn.min(p.density), mx.max(p.density))
        });
    let density_range = (max_d - min_d).max(1.0);

    for par in &res.particles {
        let px = to_px(par.position.x as f32, par.position.z as f32);
        let t = (par.speed as f32 / max_speed).clamp(0.0, 1.0);
        // Blue (slow) → cyan → green → yellow → red (fast): HSV-ish gradient.
        let color = speed_color(t);
        // Dot radius: scale with density (denser = slightly larger).
        let density_t = ((par.density - min_d) / density_range) as f32;
        let radius = 2.0 + density_t * 3.0;
        painter.circle_filled(px, radius, color);
    }

    // Box border.
    painter.rect_stroke(
        egui::Rect::from_min_max(to_px(0.0, 0.0), to_px(box_s, box_s)),
        0.0,
        egui::Stroke::new(1.0, egui::Color32::from_gray(80)),
    );

    // Stats grid below the painter.
    ui.add_space(4.0);
    egui::Grid::new("fluids_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(ui, "particles", format!("{}", res.particles.len()));
            row(ui, "steps executed", format!("{}", res.steps_done));
            row(
                ui,
                "mean density (kg/m³)",
                format!("{:.2}", res.mean_density),
            );
            row(ui, "max speed (m/s)", format!("{max_speed:.4}"));
        });
}

/// Map a normalised speed `t ∈ [0, 1]` to a colour.
///
/// Blue (t=0) → cyan → green → yellow → red (t=1), a simple HSV warm ramp.
fn speed_color(t: f32) -> egui::Color32 {
    // Four stops: blue / cyan / yellow / red.
    let stops: [(f32, [u8; 3]); 4] = [
        (0.0, [40, 80, 220]),
        (0.33, [40, 200, 200]),
        (0.66, [230, 200, 40]),
        (1.0, [220, 40, 40]),
    ];
    // Find the two surrounding stops and lerp.
    for w in stops.windows(2) {
        let (t0, c0) = w[0];
        let (t1, c1) = w[1];
        if t <= t1 {
            let f = if (t1 - t0).abs() < 1e-6 {
                0.0
            } else {
                (t - t0) / (t1 - t0)
            };
            let lerp = |a: u8, b: u8| -> u8 {
                (a as f32 + f * (b as f32 - a as f32)).clamp(0.0, 255.0) as u8
            };
            return egui::Color32::from_rgb(
                lerp(c0[0], c1[0]),
                lerp(c0[1], c1[1]),
                lerp(c0[2], c1[2]),
            );
        }
    }
    egui::Color32::from_rgb(220, 40, 40)
}

// ---------------------------------------------------------------------------
// Tests (headless_ui_tests + unit tests, mirroring sensors_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        let mut s = FluidsWorkbenchState::default();
        s.agent_set("viscosity μ (Pa·s)", &AgentValue::Float(0.5))
            .unwrap();
        assert_eq!(s.params.viscosity, 0.5);
        s.agent_set("number of steps", &AgentValue::Int(50))
            .unwrap();
        assert_eq!(s.params.num_steps, 50);
        // Unknown caption -> Err (no panic).
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into a numeric field) -> Err.
        assert!(s
            .agent_set("viscosity μ (Pa·s)", &AgentValue::Str("thick".into()))
            .is_err());
        // Out-of-range (steps > 500) -> Err, field untouched.
        assert!(s
            .agent_set("number of steps", &AgentValue::Int(9999))
            .is_err());
        assert_eq!(
            s.params.num_steps, 50,
            "rejected set leaves field untouched"
        );
    }

    #[test]
    fn default_run_succeeds_and_particles_in_box() {
        let s = FluidsWorkbenchState::default();
        let res = s.run().expect("default SPH run should succeed");
        assert!(!res.particles.is_empty(), "should have > 0 particles");
        let box_s = s.params.box_size;
        for par in &res.particles {
            assert!(
                par.position.x >= -1e-9 && par.position.x <= box_s + 1e-9,
                "particle X out of box: {}",
                par.position.x
            );
            assert!(
                par.position.z >= -1e-9 && par.position.z <= box_s + 1e-9,
                "particle Z out of box: {}",
                par.position.z
            );
        }
    }

    #[test]
    fn default_run_particle_count_matches_n_cubed() {
        let s = FluidsWorkbenchState::default();
        let res = s.run().expect("default run should succeed");
        let expected = s.params.n_per_axis.pow(3);
        assert_eq!(
            res.particles.len(),
            expected,
            "particle count must equal n_per_axis³"
        );
    }

    #[test]
    fn mean_density_is_positive_finite() {
        let s = FluidsWorkbenchState::default();
        let res = s.run().expect("run should succeed");
        assert!(
            res.mean_density.is_finite() && res.mean_density > 0.0,
            "mean density must be finite and > 0, got {}",
            res.mean_density
        );
    }

    #[test]
    fn steps_done_matches_param() {
        let s = FluidsWorkbenchState::default();
        let res = s.run().expect("run should succeed");
        assert_eq!(
            res.steps_done, s.params.num_steps,
            "steps_done must equal params.num_steps"
        );
    }

    // Degenerate-param tests — must return Err, NOT panic.

    #[test]
    fn zero_smoothing_length_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.smoothing_length = 0.0;
        assert!(s.run().is_err(), "h = 0 must return Err, not panic");
    }

    #[test]
    fn negative_smoothing_length_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.smoothing_length = -0.1;
        assert!(s.run().is_err(), "h < 0 must return Err, not panic");
    }

    #[test]
    fn nan_smoothing_length_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.smoothing_length = f64::NAN;
        assert!(s.run().is_err(), "h = NaN must return Err, not panic");
    }

    #[test]
    fn zero_dt_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.dt = 0.0;
        assert!(s.run().is_err(), "dt = 0 must return Err, not panic");
    }

    #[test]
    fn zero_n_per_axis_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.n_per_axis = 0;
        assert!(
            s.run().is_err(),
            "n_per_axis = 0 must return Err, not panic"
        );
    }

    #[test]
    fn zero_box_size_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.box_size = 0.0;
        assert!(s.run().is_err(), "box_size = 0 must return Err, not panic");
    }

    #[test]
    fn zero_num_steps_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.num_steps = 0;
        assert!(s.run().is_err(), "num_steps = 0 must return Err, not panic");
    }

    #[test]
    fn negative_viscosity_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.viscosity = -1.0;
        assert!(s.run().is_err(), "viscosity < 0 must return Err, not panic");
    }

    #[test]
    fn non_positive_rest_density_returns_err() {
        let mut s = FluidsWorkbenchState::default();
        s.params.rest_density = 0.0;
        assert!(
            s.run().is_err(),
            "rest_density = 0 must return Err, not panic"
        );
        s.params.rest_density = -500.0;
        assert!(
            s.run().is_err(),
            "rest_density < 0 must return Err, not panic"
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
            draw_fluids_workbench(app, ctx);
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
        assert!(!app.show_fluids_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_fluids_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fluids_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fluids_workbench = true;
        let res = app.fluids.run().expect("run should succeed");
        app.fluids.result = Some(res);
        app.fluids.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_fluids_workbench = true;
        // Trigger an error state.
        app.fluids.params.smoothing_length = 0.0;
        let result = app.fluids.run();
        app.fluids.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.fluids.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_fluids_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();

        // 5 SPH params + 2 integration + 2 scene = 9 DragValues minimum.
        assert!(
            spin_buttons.len() >= 9,
            "expected at least 9 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check specific captions are present as named accessibility nodes.
        for caption in [
            "smoothing length h (m)",
            "viscosity \u{03BC} (Pa\u{00B7}s)",
            "gravity |g| (m/s\u{00B2})",
            "rest density \u{03C1}\u{2080} (kg/m\u{00B3})",
            "time step dt (s)",
            "number of steps",
            "particles per axis N",
            "box size (m)",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }

        // The Run button must be named.
        assert!(
            nodes.iter().any(|(_, n)| {
                n.role() == Role::Button && n.name().is_some_and(|s| s.contains("Run simulation"))
            }),
            "the Run simulation button must be a named, invokable node"
        );
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption (egui clears a DragValue's own Name), so an AI / screen reader
        // can find the control by its caption text. Beyond merely being
        // non-empty, each `labelled_by` target must RESOLVE to the caption node
        // — i.e. the spin button is correctly associated with a real named
        // label, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_fluids_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        // 5 SPH params + 2 integration + 2 scene = 9 DragValues, all
        // unconditionally rendered (no mode gating).
        assert!(
            spin_buttons.len() >= 9,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by its caption (AI-drivable name)"
        );
        // Each spin button's labelled_by target resolves to a named caption node.
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );
        // A couple of captions must exist as named nodes in the a11y tree.
        for caption in ["smoothing length h (m)", "number of steps", "box size (m)"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn degenerate_h_shows_error_not_panic() {
        // When h ≤ 0 the workbench must surface the error in-panel, not panic.
        let mut state = FluidsWorkbenchState::default();
        state.params.smoothing_length = 0.0;
        assert!(state.run().is_err(), "h = 0 must produce Err, not panic");
        state.params.smoothing_length = -0.05;
        assert!(state.run().is_err(), "h < 0 must produce Err, not panic");
    }

    #[test]
    fn agent_bridge_fluids_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "fluids" }`:
        //   1. TabKind::from_id("fluids") → Some(TabKind::Fluids)
        //   2. set_workbench_flag(app, "fluids", true) → show_fluids_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup.
        assert_eq!(
            TabKind::from_id("fluids"),
            Some(TabKind::Fluids),
            "\"fluids\" must resolve to TabKind::Fluids"
        );
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("FLUIDS"), Some(TabKind::Fluids));
        assert_eq!(TabKind::from_id("  fluids  "), Some(TabKind::Fluids));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_fluids_workbench);
        set_workbench_flag(&mut app, "fluids", true);
        assert!(
            app.show_fluids_workbench,
            "set_workbench_flag(\"fluids\", true) must set show_fluids_workbench"
        );
        set_workbench_flag(&mut app, "fluids", false);
        assert!(!app.show_fluids_workbench);
    }
}
