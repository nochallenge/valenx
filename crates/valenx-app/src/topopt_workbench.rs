//! The right-side **Topology Optimization** workbench — a native, in-house
//! **SIMP** (Solid Isotropic Material with Penalization) generative
//! structural-design sandbox. The user picks a rectangular design domain
//! (`nx × ny`), a load case (the classic **MBB beam** or a **cantilever**), a
//! target **volume fraction** and the SIMP knobs (penalty `p`, filter radius
//! `r_min`, max iterations); pressing **Run optimization** drives the SIMP loop
//! to a **minimum-compliance** structure — the iconic organic truss generative
//! design produces.
//!
//! It is a front-end over the pure [`valenx_topopt`] crate, which does the
//! physics: bilinear Q4 plane-stress elements, an Optimality-Criteria density
//! update, a sensitivity filter, and the FE solves on the valenx-fem **faer**
//! sparse-Cholesky backend. This workbench owns only the GUI: the controls, the
//! **density heat-map** render (black = solid, white = void), the iteration
//! **animation** (the snapshot fields are replayed frame by frame so the
//! structure is seen forming), and the compliance/volume read-outs.
//!
//! Mirrors the other real-time workbenches (`morphogenesis_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_topopt_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"topopt"`. The bridge
//! can set the controls (`agent_set` / `agent_control_names`), read a status
//! line (`agent_readout`), and fire the optimisation via the RunCommand id
//! `topopt.run`.

use eframe::egui;

use valenx_topopt::{run_with_snapshots, LoadCase, TopOptParams, TopOptResult};

use crate::ValenxApp;

/// Frames each iteration snapshot is held on screen while animating (so the
/// growth is visible rather than instant). Playback advances one snapshot per
/// this-many repaints.
const FRAMES_PER_SNAPSHOT: u32 = 2;

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the Topology Optimization workbench: the input
/// parameters, the latest result (if any) and the animation cursor.
pub struct TopOptWorkbenchState {
    /// Grid width (elements along x).
    pub nx: u32,
    /// Grid height (elements along y).
    pub ny: u32,
    /// Target volume fraction `∈ (0, 1]` (the material budget).
    pub vol_frac: f64,
    /// SIMP penalty exponent `p`.
    pub penalty: f64,
    /// Sensitivity-filter radius `r_min` (element units).
    pub r_min: f64,
    /// Maximum optimisation iterations.
    pub max_iter: u32,
    /// Selected load case preset.
    pub load_case: LoadCase,
    /// The most recent optimisation result (the density field + compliance
    /// history + per-iteration snapshots). `None` until the first run.
    pub result: Option<TopOptResult>,
    /// Animation cursor: index into the result's snapshots currently shown.
    /// Advances toward the last snapshot while `playing`.
    pub anim_idx: usize,
    /// Frame sub-counter so each snapshot is held for [`FRAMES_PER_SNAPSHOT`].
    anim_frame: u32,
    /// Whether the iteration animation is currently playing.
    pub playing: bool,
}

impl Default for TopOptWorkbenchState {
    fn default() -> Self {
        let p = TopOptParams::default();
        Self {
            nx: p.nx as u32,
            ny: p.ny as u32,
            vol_frac: p.vol_frac,
            penalty: p.penalty,
            r_min: p.r_min,
            max_iter: p.max_iter as u32,
            load_case: p.load_case,
            result: None,
            anim_idx: 0,
            anim_frame: 0,
            playing: false,
        }
    }
}

impl TopOptWorkbenchState {
    /// Assemble the current controls into a clamped [`TopOptParams`].
    fn params(&self) -> TopOptParams {
        TopOptParams::new(
            self.nx as usize,
            self.ny as usize,
            self.vol_frac,
            self.penalty,
            self.r_min,
            self.max_iter as usize,
            self.load_case,
        )
    }

    /// Run the SIMP optimisation now (capturing per-iteration snapshots) and
    /// arm the animation from the first snapshot. Factored out so the in-panel
    /// **Run optimization** button and the `topopt.run` bridge id share one
    /// path.
    fn run_now(&mut self) {
        let res = run_with_snapshots(&self.params());
        self.anim_idx = 0;
        self.anim_frame = 0;
        // Animate only when there is more than one snapshot to show.
        self.playing = res.snapshots.len() > 1;
        self.result = Some(res);
    }

    /// The density field currently being displayed: the snapshot at `anim_idx`
    /// while animating, else the final field. `None` before the first run.
    fn display_density(&self) -> Option<(usize, usize, &[f64])> {
        let r = self.result.as_ref()?;
        if !r.snapshots.is_empty() {
            let i = self.anim_idx.min(r.snapshots.len() - 1);
            Some((r.nx, r.ny, &r.snapshots[i]))
        } else {
            Some((r.nx, r.ny, &r.density))
        }
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (returned by `ListControls`). Order follows the form.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "Grid nx",
            "Grid ny",
            "Volume fraction",
            "Penalty p",
            "Filter radius",
            "Max iterations",
            "Load case",
        ]
    }

    /// Set one labelled control by its caption, for the agent `SetControl`
    /// bridge. Fail-loud on an unknown caption / wrong type / out-of-range; no
    /// state is written on error and nothing panics. Numeric controls read
    /// [`crate::agent_commands::AgentValue`] as `f64`/`i64`; `Load case` reads a
    /// string and must name a known preset (`mbb` / `cantilever`).
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        match name {
            "Grid nx" => {
                let n = value.as_i64()?;
                if !(4..=200).contains(&n) {
                    return Err(format!("Grid nx must be in 4..=200, got {n}"));
                }
                self.nx = n as u32;
            }
            "Grid ny" => {
                let n = value.as_i64()?;
                if !(4..=200).contains(&n) {
                    return Err(format!("Grid ny must be in 4..=200, got {n}"));
                }
                self.ny = n as u32;
            }
            "Volume fraction" => {
                let v = value.as_f64()?;
                if !(0.05..=1.0).contains(&v) {
                    return Err(format!("Volume fraction must be in 0.05..=1.0, got {v}"));
                }
                self.vol_frac = v;
            }
            "Penalty p" => {
                let v = value.as_f64()?;
                if !(1.0..=8.0).contains(&v) {
                    return Err(format!("Penalty p must be in 1.0..=8.0, got {v}"));
                }
                self.penalty = v;
            }
            "Filter radius" => {
                let v = value.as_f64()?;
                if !(1.0..=16.0).contains(&v) {
                    return Err(format!("Filter radius must be in 1.0..=16.0, got {v}"));
                }
                self.r_min = v;
            }
            "Max iterations" => {
                let n = value.as_i64()?;
                if !(1..=1000).contains(&n) {
                    return Err(format!("Max iterations must be in 1..=1000, got {n}"));
                }
                self.max_iter = n as u32;
            }
            "Load case" => {
                let s = value.as_str()?;
                match LoadCase::from_id(s) {
                    Some(lc) => self.load_case = lc,
                    None => {
                        return Err(format!(
                            "Load case must be one of mbb/cantilever, got {s:?}"
                        ))
                    }
                }
            }
            other => return Err(format!("unknown topopt control: {other:?}")),
        }
        Ok(())
    }

    /// The current readout text for the agent `ReadReadout` bridge: the inputs
    /// plus the optimisation outcome (iterations, final compliance, achieved
    /// volume) once a run exists. `Some` once optimised, `None` before the
    /// first run (the "not run yet" surface).
    pub fn agent_readout(&self) -> Option<String> {
        let r = self.result.as_ref()?;
        Some(format!(
            "Topology optimization \u{00B7} {} \u{00B7} {}x{} grid \u{00B7} vol target {:.2} \
             \u{00B7} p={:.1} \u{00B7} r_min={:.1} \u{00B7} {} iters \u{00B7} compliance {:.4} \
             \u{00B7} achieved vol {:.3}",
            self.load_case.label(),
            r.nx,
            r.ny,
            self.vol_frac,
            self.penalty,
            self.r_min,
            r.iterations,
            r.final_compliance(),
            r.volume_fraction,
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (run optimisation)
// ---------------------------------------------------------------------------

/// Run the optimisation (the in-panel **Run optimization** action). Factored
/// out so the button and the `topopt.run` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.topopt.run_now();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Topology Optimization workbench. A no-op unless toggled on via
/// View → Topology Optimization.
pub fn draw_topopt_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_topopt_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_topopt_workbench",
        "Topology Optimization (SIMP generative design)",
        topopt_workbench_body,
    );
    if close {
        app.show_topopt_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn topopt_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| topopt_workbench_body_inner(app, ui));
}

fn topopt_workbench_body_inner(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house SIMP topology optimization \u{00B7} generative structural design [specify a \
             rectangular domain, a load case and a volume budget \u{2014} the solver iteratively \
             strips material to the stiffest (minimum-compliance) truss; FE solves run on the \
             valenx-fem faer sparse backend]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let s = &mut app.topopt;

    // --- Load-case preset selector -----------------------------------------
    ui.horizontal(|ui| {
        let lbl = ui.label("Load case");
        let mut chosen: Option<LoadCase> = None;
        egui::ComboBox::from_id_source("topopt_load_case")
            .selected_text(s.load_case.label())
            .show_ui(ui, |ui| {
                for lc in [LoadCase::MbbBeam, LoadCase::Cantilever] {
                    if ui
                        .selectable_label(s.load_case == lc, lc.label())
                        .on_hover_text(load_case_hint(lc))
                        .clicked()
                    {
                        chosen = Some(lc);
                    }
                }
            })
            .response
            .labelled_by(lbl.id);
        if let Some(lc) = chosen {
            s.load_case = lc;
        }
    });

    // --- Numeric controls ---------------------------------------------------
    egui::Grid::new("topopt_controls")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            // Grid nx.
            let lbl = ui.label("Grid nx");
            ui.add(egui::DragValue::new(&mut s.nx).speed(1).range(4..=200))
                .labelled_by(lbl.id)
                .on_hover_text("Number of elements across the domain width (x).");
            ui.end_row();

            // Grid ny.
            let lbl = ui.label("Grid ny");
            ui.add(egui::DragValue::new(&mut s.ny).speed(1).range(4..=200))
                .labelled_by(lbl.id)
                .on_hover_text("Number of elements across the domain height (y).");
            ui.end_row();

            // Volume fraction.
            let lbl = ui.label("Volume fraction");
            ui.add(
                egui::DragValue::new(&mut s.vol_frac)
                    .speed(0.01)
                    .range(0.05..=1.0)
                    .max_decimals(3),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Target material budget: the fraction of the domain left solid.");
            ui.end_row();

            // Penalty p.
            let lbl = ui.label("Penalty p");
            ui.add(
                egui::DragValue::new(&mut s.penalty)
                    .speed(0.1)
                    .range(1.0..=8.0)
                    .max_decimals(2),
            )
            .labelled_by(lbl.id)
            .on_hover_text("SIMP penalisation exponent; 3 is standard (pushes densities to 0/1).");
            ui.end_row();

            // Filter radius.
            let lbl = ui.label("Filter radius");
            ui.add(
                egui::DragValue::new(&mut s.r_min)
                    .speed(0.1)
                    .range(1.0..=16.0)
                    .max_decimals(2),
            )
            .labelled_by(lbl.id)
            .on_hover_text(
                "Sensitivity-filter radius (elements); larger = thicker members, no checkerboard.",
            );
            ui.end_row();

            // Max iterations.
            let lbl = ui.label("Max iterations");
            ui.add(
                egui::DragValue::new(&mut s.max_iter)
                    .speed(1.0)
                    .range(1..=1000),
            )
            .labelled_by(lbl.id)
            .on_hover_text("Upper bound on SIMP iterations (stops early once converged).");
            ui.end_row();
        });

    // --- Run / replay -------------------------------------------------------
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        if ui
            .button("\u{25B6} Run optimization")
            .on_hover_text(
                "Run the SIMP loop to a minimum-compliance structure and animate it forming.",
            )
            .clicked()
        {
            s.run_now();
        }
        // Replay the iteration animation from the start (if snapshots exist).
        let can_replay = s.result.as_ref().is_some_and(|r| r.snapshots.len() > 1);
        if ui
            .add_enabled(can_replay, egui::Button::new("\u{21BA} Replay"))
            .on_hover_text("Replay the iteration animation from the first step.")
            .clicked()
        {
            s.anim_idx = 0;
            s.anim_frame = 0;
            s.playing = true;
        }
    });

    // --- Drive + render the animation ---------------------------------------
    advance_animation(s, ui.ctx());

    ui.add_space(6.0);
    ui.separator();
    draw_result(s, ui);
}

/// One-line description of each load case (combo tooltips).
fn load_case_hint(lc: LoadCase) -> &'static str {
    match lc {
        LoadCase::MbbBeam => {
            "MBB beam (half-model): top-left point load, symmetry on the left, roller bottom-right."
        }
        LoadCase::Cantilever => {
            "Cantilever: left edge fully clamped, point load at the right-edge middle."
        }
    }
}

/// Advance the iteration-replay cursor while `playing`, requesting a repaint so
/// the structure visibly forms. Stops (and parks on the final field) once the
/// last snapshot is reached.
fn advance_animation(s: &mut TopOptWorkbenchState, ctx: &egui::Context) {
    if !s.playing {
        return;
    }
    let Some(last) = s
        .result
        .as_ref()
        .map(|r| r.snapshots.len().saturating_sub(1))
    else {
        s.playing = false;
        return;
    };
    if s.anim_idx >= last {
        s.anim_idx = last;
        s.playing = false;
        return;
    }
    s.anim_frame += 1;
    if s.anim_frame >= FRAMES_PER_SNAPSHOT {
        s.anim_frame = 0;
        s.anim_idx += 1;
    }
    ctx.request_repaint();
}

// ---------------------------------------------------------------------------
// Result render: density heat-map + compliance sparkline + readout
// ---------------------------------------------------------------------------

fn draw_result(s: &TopOptWorkbenchState, ui: &mut egui::Ui) {
    let Some((nx, ny, density)) = s.display_density() else {
        ui.label(
            egui::RichText::new(
                "No result yet \u{2014} set the domain, load case and volume budget, then press \
                 Run optimization.",
            )
            .italics()
            .weak(),
        );
        return;
    };

    // Status line: iteration progress (while animating) + outcome.
    if let Some(r) = s.result.as_ref() {
        let shown = if r.snapshots.is_empty() {
            r.iterations
        } else {
            s.anim_idx.min(r.snapshots.len().saturating_sub(1)) + 1
        };
        ui.label(
            egui::RichText::new(format!(
                "Iteration {shown}/{}   \u{00B7}   compliance {:.4}   \u{00B7}   achieved volume \
                 {:.3} (target {:.2})",
                r.iterations,
                r.final_compliance(),
                r.volume_fraction,
                s.vol_frac,
            ))
            .strong()
            .color(egui::Color32::from_rgb(150, 200, 230)),
        );
    }

    ui.add_space(4.0);
    ui.label(egui::RichText::new("Optimized density (black = solid, white = void)").strong());
    draw_density_heatmap(nx, ny, density, ui);

    // Compliance-history sparkline.
    if let Some(r) = s.result.as_ref() {
        if r.compliance_history.len() > 1 {
            ui.add_space(8.0);
            ui.label(egui::RichText::new("Compliance history (objective \u{2193})").strong());
            draw_compliance_plot(&r.compliance_history, ui);
        }
    }
}

/// Render the density field as a heat-map: each element a filled cell whose
/// grey level goes from white (void, `x≈0`) to black (solid, `x≈1`). Row `y=0`
/// is the bottom of the domain, so the grid is drawn flipped vertically (screen
/// y grows downward).
fn draw_density_heatmap(nx: usize, ny: usize, density: &[f64], ui: &mut egui::Ui) {
    let available = ui.available_size();
    let w = available.x.clamp(220.0, 680.0);
    // Keep the element aspect square: height follows the grid ratio.
    let cell = w / nx.max(1) as f32;
    let h = (cell * ny as f32).clamp(60.0, 460.0);
    let (rect, _) = ui.allocate_exact_size(egui::vec2(w, h), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    // Background = a mid frame so a fully-void result is still visible.
    painter.rect_filled(rect, 2.0, egui::Color32::from_rgb(235, 237, 240));

    if nx == 0 || ny == 0 || density.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "empty",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let cw = rect.width() / nx as f32;
    let ch = rect.height() / ny as f32;
    for ey in 0..ny {
        for ex in 0..nx {
            let d = density[ey * nx + ex].clamp(0.0, 1.0) as f32;
            // Solid (d=1) → black (0), void (d=0) → white (255).
            let g = ((1.0 - d) * 255.0) as u8;
            // Flip y so element row 0 is at the bottom of the panel.
            let sx = rect.left() + ex as f32 * cw;
            let sy = rect.top() + (ny - 1 - ey) as f32 * ch;
            let cell_rect =
                egui::Rect::from_min_size(egui::pos2(sx, sy), egui::vec2(cw + 0.6, ch + 0.6));
            painter.rect_filled(cell_rect, 0.0, egui::Color32::from_gray(g));
        }
    }

    painter.rect_stroke(
        rect,
        2.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(120, 130, 145)),
    );
    painter.text(
        egui::pos2(rect.right() - 4.0, rect.bottom() - 6.0),
        egui::Align2::RIGHT_BOTTOM,
        format!("{nx}x{ny} elements"),
        egui::FontId::monospace(9.0),
        egui::Color32::from_rgb(90, 100, 115),
    );
}

/// Draw the compliance history as a simple line plot (objective vs iteration),
/// normalised to its own min/max, so the user sees the "decrease then plateau".
fn draw_compliance_plot(history: &[f64], ui: &mut egui::Ui) {
    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.clamp(220.0, 680.0), 90.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, egui::Color32::from_rgb(16, 18, 24));

    let lo = history.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = history.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let span = (hi - lo).max(1e-12);
    let inner = rect.shrink(6.0);
    let n = history.len();
    let pts: Vec<egui::Pos2> = history
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            let t = if n > 1 {
                i as f32 / (n - 1) as f32
            } else {
                0.0
            };
            let yn = ((c - lo) / span) as f32; // 0 = min (best), 1 = max
            egui::pos2(
                inner.left() + t * inner.width(),
                inner.top() + yn * inner.height(),
            )
        })
        .collect();
    painter.add(egui::Shape::line(
        pts,
        egui::Stroke::new(1.6, egui::Color32::from_rgb(120, 200, 255)),
    ));
    painter.rect_stroke(
        rect,
        2.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(50, 60, 80)),
    );
    painter.text(
        egui::pos2(rect.left() + 4.0, rect.top() + 3.0),
        egui::Align2::LEFT_TOP,
        format!("c: {hi:.3} \u{2192} {lo:.3}"),
        egui::FontId::monospace(9.0),
        egui::Color32::from_gray(170),
    );
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_state_matches_param_defaults() {
        let s = TopOptWorkbenchState::default();
        let p = TopOptParams::default();
        assert_eq!(s.nx as usize, p.nx);
        assert_eq!(s.ny as usize, p.ny);
        assert_eq!(s.load_case, LoadCase::MbbBeam);
        assert!(s.result.is_none());
        assert!(!s.playing);
    }

    #[test]
    fn agent_set_numeric_controls() {
        use crate::agent_commands::AgentValue;
        let mut s = TopOptWorkbenchState::default();
        s.agent_set("Grid nx", &AgentValue::Int(80)).expect("nx");
        assert_eq!(s.nx, 80);
        s.agent_set("Grid ny", &AgentValue::Int(40)).expect("ny");
        assert_eq!(s.ny, 40);
        s.agent_set("Volume fraction", &AgentValue::Float(0.35))
            .expect("vf");
        assert!((s.vol_frac - 0.35).abs() < 1e-9);
        s.agent_set("Penalty p", &AgentValue::Float(4.0))
            .expect("p");
        assert!((s.penalty - 4.0).abs() < 1e-9);
        s.agent_set("Filter radius", &AgentValue::Float(2.5))
            .expect("r");
        assert!((s.r_min - 2.5).abs() < 1e-9);
        s.agent_set("Max iterations", &AgentValue::Int(40))
            .expect("it");
        assert_eq!(s.max_iter, 40);
    }

    #[test]
    fn agent_set_load_case() {
        use crate::agent_commands::AgentValue;
        let mut s = TopOptWorkbenchState::default();
        s.agent_set("Load case", &AgentValue::Str("cantilever".to_string()))
            .expect("valid load case");
        assert_eq!(s.load_case, LoadCase::Cantilever);
    }

    #[test]
    fn agent_set_rejects_out_of_range_and_unknown() {
        use crate::agent_commands::AgentValue;
        let mut s = TopOptWorkbenchState::default();
        assert!(s.agent_set("Grid nx", &AgentValue::Int(1)).is_err());
        assert!(s.agent_set("Grid nx", &AgentValue::Int(9999)).is_err());
        assert!(s
            .agent_set("Volume fraction", &AgentValue::Float(2.0))
            .is_err());
        assert!(s.agent_set("Penalty p", &AgentValue::Float(0.1)).is_err());
        assert!(s
            .agent_set("Filter radius", &AgentValue::Float(0.0))
            .is_err());
        assert!(s
            .agent_set("Load case", &AgentValue::Str("nope".to_string()))
            .is_err());
        assert!(s.agent_set("bogus", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn readout_is_none_until_run_then_reports_outcome() {
        let mut s = TopOptWorkbenchState::default();
        assert!(s.agent_readout().is_none(), "no readout before a run");
        // Use a small fast problem for the test.
        s.nx = 30;
        s.ny = 15;
        s.max_iter = 12;
        s.run_now();
        let r = s.agent_readout().expect("readout after run");
        assert!(r.contains("Topology optimization"), "readout: {r}");
        assert!(r.contains("compliance"), "readout names compliance: {r}");
        assert!(r.contains("achieved vol"), "readout names volume: {r}");
    }

    #[test]
    fn run_now_arms_animation_and_produces_a_field() {
        let mut s = TopOptWorkbenchState {
            nx: 30,
            ny: 15,
            max_iter: 12,
            ..Default::default()
        };
        s.run_now();
        let r = s.result.as_ref().expect("a result after run");
        assert!(!r.density.is_empty());
        assert!(!r.snapshots.is_empty(), "snapshots captured for animation");
        assert_eq!(s.anim_idx, 0, "animation starts at the first snapshot");
        // display_density should yield the first snapshot while playing.
        let (nx, ny, _) = s.display_density().expect("display field");
        assert_eq!((nx, ny), (30, 15));
    }

    #[test]
    fn run_bridge_helper_runs_through_app() {
        let mut app = ValenxApp::default();
        app.topopt.nx = 28;
        app.topopt.ny = 14;
        app.topopt.max_iter = 10;
        run(&mut app);
        assert!(
            app.topopt.result.is_some(),
            "the topopt.run bridge helper should produce a result"
        );
    }

    #[test]
    fn animation_advances_then_parks_on_final() {
        let ctx = egui::Context::default();
        let mut s = TopOptWorkbenchState {
            nx: 24,
            ny: 12,
            max_iter: 10,
            ..Default::default()
        };
        s.run_now();
        let last = s.result.as_ref().unwrap().snapshots.len() - 1;
        // Pump frames; eventually anim_idx reaches the last snapshot and stops.
        for _ in 0..(last + 2) * FRAMES_PER_SNAPSHOT as usize + 4 {
            advance_animation(&mut s, &ctx);
        }
        assert_eq!(s.anim_idx, last, "animation parks on the final snapshot");
        assert!(!s.playing, "animation stops at the end");
    }

    #[test]
    fn control_names_are_listed() {
        let names = TopOptWorkbenchState::agent_control_names();
        for c in [
            "Grid nx",
            "Grid ny",
            "Volume fraction",
            "Penalty p",
            "Filter radius",
            "Max iterations",
            "Load case",
        ] {
            assert!(names.contains(&c), "missing control name {c}");
        }
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
            draw_topopt_workbench(app, ctx);
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
        assert!(!app.show_topopt_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_topopt_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_topopt_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_topopt_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 6,
            "expected the six numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );

        for caption in [
            "Grid nx",
            "Grid ny",
            "Volume fraction",
            "Penalty p",
            "Filter radius",
            "Max iterations",
            "Load case",
        ] {
            assert!(
                has_named_node(&nodes, caption),
                "caption '{caption}' must be a named node in the a11y tree"
            );
        }
    }
}
