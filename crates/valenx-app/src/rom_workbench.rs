//! The right-side **ROM Workbench** panel — a native front-end over the
//! in-house `valenx-rom` crate (reduced-order & surrogate modelling).
//!
//! This workbench demonstrates **Proper Orthogonal Decomposition (POD)** — the
//! energy-optimal low-rank basis of a snapshot matrix via its SVD (see
//! [`valenx_rom::PodBasis`]). The user dials a **synthetic 1-D parametric
//! field** `f(x; μ)` — a sum of moving Gaussian bumps whose centres sweep with a
//! varying parameter μ — sampled on a spatial grid into a snapshot matrix
//! (columns = parameter samples). Clicking **Run** fits the POD basis, reads off
//! the singular-value energy spectrum, truncates to a chosen rank `k`, and
//! reconstructs the field, reporting the captured-energy fraction and the
//! reconstruction RMSE.
//!
//! Mirrors the other workbenches (`ocean_workbench`, `fluids_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_rom_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"rom"` (see
//! [`crate::project_tabs::TabKind`] / [`crate::agent_commands`]).
//!
//! Two painter views are drawn: (a) the singular-value / cumulative-energy decay
//! versus mode index (log-y), and (b) a chosen snapshot column overlaid with its
//! `k`-mode POD reconstruction.
//!
//! Honesty: `valenx-rom` is **research / educational-grade textbook ROM** — POD
//! here is an exact SVD construction (validated to machine precision against
//! analytic fields in-crate), but the *field generator* in this panel is a
//! synthetic teaching example, not data from a physical solve. `valenx-rom` also
//! provides DMD (dynamic mode decomposition) and DEIM (hyper-reduction); this
//! first-cut workbench surfaces POD. Every error from `valenx-rom` surfaces
//! verbatim — the workbench never invents a number, and degenerate parameters
//! (e.g. `resolution = 0`, `rank = 0`, `rank > #snapshots`) show an in-panel
//! error, NOT a panic.

use eframe::egui;
use valenx_rom::{PodBasis, Snapshots};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Editable synthetic-field + POD parameters shown in the workbench.
#[derive(Clone, Debug)]
pub struct RomParams {
    // --- Synthetic 1-D parametric field f(x; μ) ---
    /// Spatial resolution `n` — number of grid points the field is sampled on
    /// (the snapshot state dimension). Must be ≥ 1.
    pub resolution: usize,
    /// Number of snapshots `m` — parameter samples `μ_0 … μ_{m-1}` swept across
    /// the field (the snapshot column count). Must be ≥ 1.
    pub num_snapshots: usize,
    /// Number of Gaussian bumps summed into each snapshot — must be ≥ 1.
    pub num_bumps: usize,
    /// Feature width `σ` (fraction of the domain) of each Gaussian bump — must
    /// be > 0. Narrower bumps need more modes to resolve.
    pub feature_width: f64,
    /// How far (fraction of the domain) the bump centres travel as μ sweeps from
    /// its first to its last sample. Larger travel → a richer (higher-rank)
    /// snapshot set.
    pub travel: f64,
    /// Additive Gaussian-like noise amplitude (deterministic pseudo-noise so the
    /// run is reproducible), ≥ 0. Noise lifts the singular-value tail.
    pub noise: f64,

    // --- POD truncation ---
    /// Truncation rank `k` — how many POD modes the reconstruction keeps. Must
    /// be ≥ 1 and ≤ the number of significant singular values.
    pub rank: usize,
    /// Which snapshot column to overlay against its reconstruction in view (b).
    pub view_snapshot: usize,
}

impl Default for RomParams {
    fn default() -> Self {
        Self {
            resolution: 80,
            num_snapshots: 40,
            num_bumps: 2,
            feature_width: 0.08,
            travel: 0.6,
            noise: 0.0,
            rank: 4,
            view_snapshot: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// Cached POD output for the painter + readouts.
#[derive(Default, Clone)]
pub struct RomResult {
    /// **All** numerically significant singular values of the snapshot matrix,
    /// descending (for the energy-spectrum painter). Length = data rank.
    pub singular_values: Vec<f64>,
    /// Cumulative captured-energy fraction `Σ_{i≤j} σ_i² / Σ σ_i²` per index `j`
    /// (length-aligned with [`Self::singular_values`]), in `[0, 1]`.
    pub cumulative_energy: Vec<f64>,
    /// The retained truncation rank `k` actually used (clamped into range).
    pub rank: usize,
    /// Fraction of total snapshot energy captured by the `k` retained modes.
    pub captured_energy: f64,
    /// Relative Frobenius reconstruction error `‖X − UₖUₖᵀX‖_F / ‖X‖_F`.
    pub reconstruction_error: f64,
    /// Absolute root-mean-square reconstruction error over every grid point of
    /// every snapshot (same units as the field).
    pub reconstruction_rmse: f64,
    /// Fewest modes whose cumulative energy reaches 99 %.
    pub modes_for_99: usize,
    /// The chosen snapshot column (the "truth" curve in view (b)), one value per
    /// grid point.
    pub snapshot_truth: Vec<f64>,
    /// The same column reconstructed from its `k`-mode projection.
    pub snapshot_recon: Vec<f64>,
    /// Index of the snapshot drawn in view (b).
    pub snapshot_index: usize,
    /// Grid coordinates `x ∈ [0, 1]` for the snapshot/reconstruction curves.
    pub x_grid: Vec<f64>,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the ROM workbench.
#[derive(Default)]
pub struct RomWorkbenchState {
    /// User-editable parameters.
    pub params: RomParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<RomResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl RomWorkbenchState {
    /// Generate the synthetic snapshot matrix `X` (state_dim × n_time) from the
    /// current parameters, as a [`Snapshots`].
    ///
    /// Column `j` is the field `f(x; μ_j)` sampled on the spatial grid, where
    /// `μ_j` sweeps the bump centres across the domain. Returns `Err` (shown
    /// in-panel) — never panics — when the user has entered degenerate values
    /// (`resolution == 0`, `num_snapshots == 0`, non-positive `feature_width`),
    /// delegating the final non-empty / finite check to `valenx-rom`'s fail-loud
    /// [`Snapshots::from_columns`].
    pub fn build_snapshots(&self) -> Result<Snapshots, String> {
        let p = &self.params;
        if p.resolution == 0 {
            return Err("spatial resolution must be >= 1, got 0".to_string());
        }
        if p.num_snapshots == 0 {
            return Err("number of snapshots must be >= 1, got 0".to_string());
        }
        if p.num_bumps == 0 {
            return Err("number of bumps must be >= 1, got 0".to_string());
        }
        if !(p.feature_width.is_finite() && p.feature_width > 0.0) {
            return Err(format!(
                "feature width must be finite and > 0, got {}",
                p.feature_width
            ));
        }
        if !p.travel.is_finite() {
            return Err(format!("travel must be finite, got {}", p.travel));
        }
        if !(p.noise.is_finite() && p.noise >= 0.0) {
            return Err(format!("noise must be finite and >= 0, got {}", p.noise));
        }

        let n = p.resolution;
        let m = p.num_snapshots;
        let sigma = p.feature_width;
        let two_sig2 = 2.0 * sigma * sigma;

        let mut columns: Vec<Vec<f64>> = Vec::with_capacity(m);
        for j in 0..m {
            // Parameter μ_j ∈ [0, 1] over the snapshot index.
            let mu = if m == 1 {
                0.0
            } else {
                j as f64 / (m - 1) as f64
            };
            let mut col = Vec::with_capacity(n);
            for i in 0..n {
                // Spatial coordinate x ∈ [0, 1].
                let x = if n == 1 {
                    0.5
                } else {
                    i as f64 / (n - 1) as f64
                };
                let mut v = 0.0_f64;
                for b in 0..p.num_bumps {
                    // Each bump starts spread across the domain and travels by
                    // `travel · μ` (with a per-bump amplitude + a small phase so
                    // multiple bumps are linearly independent in time).
                    let base = (b as f64 + 1.0) / (p.num_bumps as f64 + 1.0);
                    let dir = if b % 2 == 0 { 1.0 } else { -1.0 };
                    let centre = (base + dir * p.travel * mu).rem_euclid(1.0);
                    let amp = 1.0 + 0.5 * b as f64;
                    let d = x - centre;
                    v += amp * (-(d * d) / two_sig2).exp();
                }
                // Deterministic pseudo-noise: a fixed high-frequency function of
                // (i, j) so the run is reproducible (no RNG dependency) yet lifts
                // the singular-value tail like measurement noise would.
                if p.noise > 0.0 {
                    let phase = (i as f64) * 12.9898 + (j as f64) * 78.233;
                    let pseudo = (phase.sin() * 43758.5453).fract() * 2.0 - 1.0;
                    v += p.noise * pseudo;
                }
                col.push(v);
            }
            columns.push(col);
        }

        Snapshots::from_columns(&columns).map_err(|e| e.to_string())
    }

    /// Run the full POD pipeline: build the snapshot matrix, fit a full-rank POD
    /// basis for the energy spectrum, fit a rank-`k` basis for the truncated
    /// reconstruction, and report the captured energy + reconstruction errors.
    ///
    /// Every failure is returned as an `Err(String)` — no panics, no invented
    /// numbers. A `rank` of 0 or one exceeding the data's significant-mode count
    /// surfaces `valenx-rom`'s [`valenx_rom::RomError::InvalidRank`].
    pub fn run(&self) -> Result<RomResult, String> {
        let p = &self.params;
        let snaps = self.build_snapshots()?;

        // Full-rank basis (energy_tol = 1.0 keeps every significant mode) — its
        // retained singular values ARE the full significant spectrum, which we
        // use both to draw the decay and to compute the 99%-energy mode count.
        let full = PodBasis::fit(&snaps, 1.0).map_err(|e| e.to_string())?;
        let singular_values: Vec<f64> = full.singular_values().iter().copied().collect();

        // Cumulative captured-energy fraction per index.
        let total_energy: f64 = singular_values.iter().map(|s| s * s).sum();
        let mut cumulative_energy = Vec::with_capacity(singular_values.len());
        let mut modes_for_99 = singular_values.len();
        {
            let mut acc = 0.0_f64;
            let mut hit99 = false;
            for (idx, s) in singular_values.iter().enumerate() {
                acc += s * s;
                let frac = if total_energy > 0.0 {
                    acc / total_energy
                } else {
                    0.0
                };
                cumulative_energy.push(frac);
                if !hit99 && frac >= 0.99 {
                    modes_for_99 = idx + 1;
                    hit99 = true;
                }
            }
        }

        // Rank-k truncated basis for the reconstruction. `fit_rank` is fail-loud
        // on rank 0 or rank > significant-mode count — surfaced in-panel.
        let basis = PodBasis::fit_rank(&snaps, p.rank).map_err(|e| e.to_string())?;
        let rank = basis.rank();
        let captured_energy = basis.captured_energy();
        let reconstruction_error = basis
            .reconstruction_error(&snaps)
            .map_err(|e| e.to_string())?;

        // Absolute RMSE over the full reconstructed matrix.
        let x = snaps.matrix();
        let mut sq_sum = 0.0_f64;
        let count = x.nrows() * x.ncols();
        for col in 0..x.ncols() {
            let xc = x.column(col).into_owned();
            let a = basis.project(&xc).map_err(|e| e.to_string())?;
            let rc = basis.reconstruct(&a).map_err(|e| e.to_string())?;
            for row in 0..x.nrows() {
                let d = xc[row] - rc[row];
                sq_sum += d * d;
            }
        }
        let reconstruction_rmse = if count > 0 {
            (sq_sum / count as f64).sqrt()
        } else {
            0.0
        };

        // The chosen snapshot column + its k-mode reconstruction for view (b).
        let n_time = snaps.n_time();
        let snapshot_index = p.view_snapshot.min(n_time.saturating_sub(1));
        let truth_col = x.column(snapshot_index).into_owned();
        let a = basis.project(&truth_col).map_err(|e| e.to_string())?;
        let recon_col = basis.reconstruct(&a).map_err(|e| e.to_string())?;
        let snapshot_truth: Vec<f64> = truth_col.iter().copied().collect();
        let snapshot_recon: Vec<f64> = recon_col.iter().copied().collect();
        let n = snapshot_truth.len();
        let x_grid: Vec<f64> = (0..n)
            .map(|i| {
                if n <= 1 {
                    0.5
                } else {
                    i as f64 / (n - 1) as f64
                }
            })
            .collect();

        Ok(RomResult {
            singular_values,
            cumulative_energy,
            rank,
            captured_energy,
            reconstruction_error,
            reconstruction_rmse,
            modes_for_99,
            snapshot_truth,
            snapshot_recon,
            snapshot_index,
            x_grid,
        })
    }

    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by `ListControls`;
    /// each string matches exactly the caption the form draws.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "spatial resolution n",
            "number of snapshots m",
            "number of bumps",
            "feature width σ",
            "centre travel",
            "noise amplitude",
            "truncation rank k",
            "view snapshot index",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. Fail-loud: an unknown caption or a value of the
    /// wrong type / out of range returns `Err(String)` — never a panic. Ranges
    /// mirror the form's `DragValue` clamps exactly (the integer counts are all
    /// `>= 1` except the zero-based `view snapshot index`).
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
            "spatial resolution n" => {
                p.resolution = ranged_int(value, 1, 2000, "spatial resolution n")?
            }
            "number of snapshots m" => {
                p.num_snapshots = ranged_int(value, 1, 2000, "number of snapshots m")?
            }
            "number of bumps" => p.num_bumps = ranged_int(value, 1, 16, "number of bumps")?,
            "feature width σ" => {
                p.feature_width = ranged(value.as_f64()?, 0.005, 0.5, "feature width σ")?
            }
            "centre travel" => p.travel = ranged(value.as_f64()?, 0.0, 1.0, "centre travel")?,
            "noise amplitude" => p.noise = ranged(value.as_f64()?, 0.0, 2.0, "noise amplitude")?,
            "truncation rank k" => p.rank = ranged_int(value, 1, 2000, "truncation rank k")?,
            "view snapshot index" => {
                p.view_snapshot = ranged_int(value, 0, 2000, "view snapshot index")?
            }
            other => return Err(format!("unknown ROM control: {other:?}")),
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the ROM workbench. A no-op unless toggled on via View → ROM.
///
/// Mirrors [`crate::ocean_workbench::draw_ocean_workbench`].
pub fn draw_rom_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_rom_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_rom_workbench",
        "Reduced-order model (POD)",
        rom_workbench_body,
    );
    if close {
        app.show_rom_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn rom_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Proper Orthogonal Decomposition (POD) — energy-optimal SVD basis of a synthetic \
             1-D parametric field f(x; \u{03BC}) of moving Gaussian bumps · valenx-rom  \
             [research / educational — exact SVD; DMD + DEIM also in valenx-rom]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.rom;
        let p = &mut s.params;

        ui.label(egui::RichText::new("Synthetic field f(x; \u{03BC})").strong());
        egui::Grid::new("rom_field_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("spatial resolution n");
                ui.add(
                    egui::DragValue::new(&mut p.resolution)
                        .speed(1)
                        .range(1..=2000),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Grid points the field is sampled on (the snapshot state \
                     dimension). Must be >= 1.",
                );
                ui.end_row();

                let lbl = ui.label("number of snapshots m");
                ui.add(
                    egui::DragValue::new(&mut p.num_snapshots)
                        .speed(1)
                        .range(1..=2000),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Parameter samples \u{03BC}_0 \u{2026} \u{03BC}_(m-1) swept across the \
                     field (the snapshot column count). Must be >= 1.",
                );
                ui.end_row();

                let lbl = ui.label("number of bumps");
                ui.add(
                    egui::DragValue::new(&mut p.num_bumps)
                        .speed(1)
                        .range(1..=16),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "How many travelling Gaussian bumps are summed into each \
                     snapshot. More bumps = a higher-rank field.",
                );
                ui.end_row();

                let lbl = ui.label("feature width \u{03C3}");
                ui.add(
                    egui::DragValue::new(&mut p.feature_width)
                        .speed(0.005)
                        .range(0.005..=0.5),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Gaussian bump width as a fraction of the domain. Must be > 0. \
                     Narrower features need more POD modes to resolve.",
                );
                ui.end_row();

                let lbl = ui.label("centre travel");
                ui.add(
                    egui::DragValue::new(&mut p.travel)
                        .speed(0.02)
                        .range(0.0..=1.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "How far the bump centres move (fraction of the domain) as \u{03BC} \
                     sweeps. More travel = a richer, higher-rank snapshot set.",
                );
                ui.end_row();

                let lbl = ui.label("noise amplitude");
                ui.add(
                    egui::DragValue::new(&mut p.noise)
                        .speed(0.005)
                        .range(0.0..=2.0),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Deterministic pseudo-noise added to every sample (>= 0). \
                     Noise lifts the singular-value tail / raises the truncation error.",
                );
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("POD truncation").strong());
        egui::Grid::new("rom_pod_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("truncation rank k");
                ui.add(egui::DragValue::new(&mut p.rank).speed(1).range(1..=2000))
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "How many POD modes the reconstruction keeps. Must be >= 1 and \
                         <= the number of significant singular values; more modes \
                         lower the reconstruction error.",
                    );
                ui.end_row();

                let lbl = ui.label("view snapshot index");
                ui.add(
                    egui::DragValue::new(&mut p.view_snapshot)
                        .speed(1)
                        .range(0..=2000),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Which snapshot column is overlaid against its k-mode \
                     reconstruction in the lower plot.",
                );
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text("Build the snapshot matrix, fit POD, and reconstruct with k modes.")
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside borrow) -------------------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.rom;
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
    draw_rom_viz(s, ui);
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.rom;
    match s.run() {
        Ok(res) => {
            s.status = format!(
                "\u{2714} k = {} modes \u{00B7} {:.3}% energy \u{00B7} RMSE {:.3e} \u{00B7} \
                 {} modes for 99%",
                res.rank,
                100.0 * res.captured_energy,
                res.reconstruction_rmse,
                res.modes_for_99,
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
// 2-D visualisation (energy spectrum + snapshot-vs-reconstruction overlay)
// ---------------------------------------------------------------------------

fn draw_rom_viz(s: &RomWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new("press \"Run\" to fit POD and visualise the energy spectrum")
                .weak(),
        );
        return;
    };

    draw_energy_spectrum(res, ui);
    ui.add_space(8.0);
    draw_reconstruction_overlay(res, ui);

    // Readouts grid below the painters.
    ui.add_space(6.0);
    egui::Grid::new("rom_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(ui, "retained rank k", format!("{}", res.rank));
            row(
                ui,
                "captured energy at k",
                format!("{:.4} %", 100.0 * res.captured_energy),
            );
            row(
                ui,
                "reconstruction RMSE (abs)",
                format!("{:.4e}", res.reconstruction_rmse),
            );
            row(
                ui,
                "reconstruction error (rel Frobenius)",
                format!("{:.4e}", res.reconstruction_error),
            );
            row(ui, "modes for 99% energy", format!("{}", res.modes_for_99));
            row(
                ui,
                "significant modes (data rank)",
                format!("{}", res.singular_values.len()),
            );
        });
}

/// View (a): singular-value decay (log-y bars) + cumulative-energy curve.
fn draw_energy_spectrum(res: &RomResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Singular-value spectrum + cumulative energy (log-y)").strong());
    ui.label(
        egui::RichText::new(
            "amber bars = \u{03C3}_i (log scale) \u{00B7} green = cumulative captured energy \
             \u{00B7} dashed = retained rank k",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 200.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.singular_values.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no singular values",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let margin = 16.0_f32;
    let inner = rect.shrink(margin);

    // Log-y range over the singular values. Guard against zeros/negatives by
    // flooring at a small fraction of the max (the crate already drops modes at
    // the noise floor, so all retained σ are strictly positive, but be safe).
    let smax = res
        .singular_values
        .iter()
        .cloned()
        .fold(0.0_f64, f64::max)
        .max(f64::MIN_POSITIVE);
    let log_max = smax.log10();
    let log_min = (smax * 1e-6).max(f64::MIN_POSITIVE).log10();
    let log_span = (log_max - log_min).max(1e-9);

    let count = res.singular_values.len();
    let bar_w = (inner.width() / count as f32).max(1.0);

    for (i, &s) in res.singular_values.iter().enumerate() {
        let s_floored = s.max(smax * 1e-6).max(f64::MIN_POSITIVE);
        let frac = ((s_floored.log10() - log_min) / log_span).clamp(0.0, 1.0) as f32;
        let x0 = inner.left() + i as f32 * bar_w;
        let h = frac * inner.height();
        let bar = egui::Rect::from_min_max(
            egui::pos2(x0 + 0.5, inner.bottom() - h),
            egui::pos2(x0 + bar_w - 0.5, inner.bottom()),
        );
        // Highlight retained modes (index < k) brighter than the dropped tail.
        let col = if i < res.rank {
            egui::Color32::from_rgb(230, 180, 70)
        } else {
            egui::Color32::from_rgb(120, 95, 45)
        };
        painter.rect_filled(bar, 0.0, col);
    }

    // Cumulative-energy curve (linear 0..1 mapped to the full height).
    if res.cumulative_energy.len() >= 2 {
        let pts: Vec<egui::Pos2> = res
            .cumulative_energy
            .iter()
            .enumerate()
            .map(|(i, &e)| {
                let nx = (i as f32 + 0.5) / count as f32;
                let ny = 1.0 - e.clamp(0.0, 1.0) as f32;
                egui::pos2(
                    inner.left() + nx * inner.width(),
                    inner.top() + ny * inner.height(),
                )
            })
            .collect();
        painter.add(egui::Shape::line(
            pts,
            egui::Stroke::new(2.0, egui::Color32::from_rgb(110, 200, 130)),
        ));
    }

    // Dashed vertical line at the retained rank k.
    if res.rank >= 1 && res.rank <= count {
        let kx = inner.left() + (res.rank as f32 / count as f32) * inner.width();
        let mut y = inner.top();
        while y < inner.bottom() {
            let y2 = (y + 5.0).min(inner.bottom());
            painter.line_segment(
                [egui::pos2(kx, y), egui::pos2(kx, y2)],
                egui::Stroke::new(1.0, egui::Color32::from_gray(150)),
            );
            y += 9.0;
        }
    }
}

/// View (b): the chosen snapshot column overlaid with its k-mode POD
/// reconstruction.
fn draw_reconstruction_overlay(res: &RomResult, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(format!(
            "Snapshot #{} vs its {}-mode POD reconstruction",
            res.snapshot_index, res.rank
        ))
        .strong(),
    );
    ui.label(
        egui::RichText::new("blue = truth f(x; \u{03BC}) \u{00B7} orange = reconstruction")
            .weak()
            .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 200.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.snapshot_truth.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "snapshot too small to plot",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let margin = 16.0_f32;
    let inner = rect.shrink(margin);

    // Shared vertical range over both curves so the overlay is comparable.
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for v in res.snapshot_truth.iter().chain(res.snapshot_recon.iter()) {
        if v.is_finite() {
            lo = lo.min(*v);
            hi = hi.max(*v);
        }
    }
    if !(lo.is_finite() && hi.is_finite()) {
        lo = 0.0;
        hi = 1.0;
    }
    let pad = ((hi - lo) * 0.08).max(1e-6);
    let y_lo = lo - pad;
    let y_hi = hi + pad;
    let y_span = (y_hi - y_lo).max(1e-9);

    let n = res.snapshot_truth.len();
    let to_px = |i: usize, v: f64| -> egui::Pos2 {
        let nx = if n <= 1 {
            0.5
        } else {
            i as f32 / (n - 1) as f32
        };
        let ny = ((y_hi - v) / y_span).clamp(0.0, 1.0) as f32;
        egui::pos2(
            inner.left() + nx * inner.width(),
            inner.top() + ny * inner.height(),
        )
    };

    // Zero reference line if the range straddles zero.
    if y_lo < 0.0 && y_hi > 0.0 {
        let zy = to_px(0, 0.0).y;
        painter.line_segment(
            [egui::pos2(inner.left(), zy), egui::pos2(inner.right(), zy)],
            egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
        );
    }

    let truth_pts: Vec<egui::Pos2> = res
        .snapshot_truth
        .iter()
        .enumerate()
        .map(|(i, &v)| to_px(i, v))
        .collect();
    let recon_pts: Vec<egui::Pos2> = res
        .snapshot_recon
        .iter()
        .enumerate()
        .map(|(i, &v)| to_px(i, v))
        .collect();

    painter.add(egui::Shape::line(
        truth_pts,
        egui::Stroke::new(2.0, egui::Color32::from_rgb(120, 190, 240)),
    ));
    painter.add(egui::Shape::line(
        recon_pts,
        egui::Stroke::new(1.5, egui::Color32::from_rgb(240, 160, 80)),
    ));
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring ocean_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    #[test]
    fn agent_set_sets_param_unknown_and_type_mismatch_err() {
        let mut s = RomWorkbenchState::default();
        s.agent_set("truncation rank k", &AgentValue::Int(5))
            .unwrap();
        assert_eq!(s.params.rank, 5);
        s.agent_set("feature width σ", &AgentValue::Float(0.1))
            .unwrap();
        assert_eq!(s.params.feature_width, 0.1);
        // Unknown caption -> Err (no panic).
        assert!(s.agent_set("no such control", &AgentValue::Int(1)).is_err());
        // Type mismatch (string into an integer field) -> Err.
        assert!(s
            .agent_set("truncation rank k", &AgentValue::Str("five".into()))
            .is_err());
        // Out-of-range (rank 0) -> Err, field untouched.
        assert!(s
            .agent_set("truncation rank k", &AgentValue::Int(0))
            .is_err());
        assert_eq!(s.params.rank, 5, "rejected set leaves field untouched");
    }

    #[test]
    fn default_run_succeeds_and_spectrum_is_populated() {
        let s = RomWorkbenchState::default();
        let res = s.run().expect("default ROM run should succeed");
        assert!(
            !res.singular_values.is_empty(),
            "spectrum must have at least one singular value"
        );
        assert_eq!(
            res.cumulative_energy.len(),
            res.singular_values.len(),
            "cumulative energy is index-aligned with the spectrum"
        );
        assert!(
            res.captured_energy.is_finite()
                && res.captured_energy >= 0.0
                && res.captured_energy <= 1.0 + 1e-9,
            "captured energy must be a fraction, got {}",
            res.captured_energy
        );
        assert!(
            res.reconstruction_rmse.is_finite() && res.reconstruction_rmse >= 0.0,
            "RMSE must be finite and >= 0"
        );
        assert!(res.modes_for_99 >= 1 && res.modes_for_99 <= res.singular_values.len());
    }

    #[test]
    fn singular_values_are_descending() {
        let s = RomWorkbenchState::default();
        let res = s.run().expect("run");
        for w in res.singular_values.windows(2) {
            assert!(
                w[0] >= w[1] - 1e-12,
                "singular values must be descending: {} < {}",
                w[0],
                w[1]
            );
        }
    }

    #[test]
    fn cumulative_energy_reaches_one() {
        let s = RomWorkbenchState::default();
        let res = s.run().expect("run");
        let last = *res.cumulative_energy.last().expect("non-empty");
        assert!(
            (last - 1.0).abs() < 1e-9,
            "cumulative energy must end at 1.0, got {last}"
        );
    }

    #[test]
    fn more_modes_lower_or_equal_reconstruction_error() {
        // The headline POD property: increasing the truncation rank never
        // increases the reconstruction error (Eckart–Young monotonicity).
        let mut few = RomWorkbenchState::default();
        few.params.rank = 2;
        let mut many = RomWorkbenchState::default();
        many.params.rank = 8;
        let e_few = few.run().expect("few-mode run").reconstruction_error;
        let e_many = many.run().expect("many-mode run").reconstruction_error;
        assert!(
            e_many <= e_few + 1e-12,
            "more modes must not increase error: {e_many} > {e_few}"
        );
    }

    #[test]
    fn more_modes_lower_or_equal_rmse() {
        // Same monotonicity, expressed in the absolute RMSE readout.
        let mut few = RomWorkbenchState::default();
        few.params.rank = 1;
        let mut many = RomWorkbenchState::default();
        many.params.rank = 6;
        let r_few = few.run().expect("few").reconstruction_rmse;
        let r_many = many.run().expect("many").reconstruction_rmse;
        assert!(
            r_many <= r_few + 1e-12,
            "more modes must not increase RMSE: {r_many} > {r_few}"
        );
    }

    #[test]
    fn captured_energy_grows_with_rank() {
        let mut few = RomWorkbenchState::default();
        few.params.rank = 1;
        let mut many = RomWorkbenchState::default();
        many.params.rank = 5;
        let c_few = few.run().expect("few").captured_energy;
        let c_many = many.run().expect("many").captured_energy;
        assert!(
            c_many >= c_few - 1e-12,
            "more modes must not lower captured energy: {c_many} < {c_few}"
        );
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_resolution_returns_err() {
        let mut s = RomWorkbenchState::default();
        s.params.resolution = 0;
        assert!(
            s.run().is_err(),
            "resolution = 0 must return Err, not panic"
        );
    }

    #[test]
    fn zero_snapshots_returns_err() {
        let mut s = RomWorkbenchState::default();
        s.params.num_snapshots = 0;
        assert!(
            s.run().is_err(),
            "num_snapshots = 0 must return Err, not panic"
        );
    }

    #[test]
    fn zero_rank_returns_err() {
        let mut s = RomWorkbenchState::default();
        s.params.rank = 0;
        assert!(s.run().is_err(), "rank = 0 must return Err, not panic");
    }

    #[test]
    fn rank_above_data_rank_returns_err() {
        // A trivially low-rank field has few significant modes; a huge rank must
        // surface InvalidRank, not panic.
        let mut s = RomWorkbenchState::default();
        s.params.rank = 1000;
        assert!(
            s.run().is_err(),
            "rank > significant modes must return Err, not panic"
        );
    }

    #[test]
    fn zero_feature_width_returns_err() {
        let mut s = RomWorkbenchState::default();
        s.params.feature_width = 0.0;
        assert!(
            s.run().is_err(),
            "feature width = 0 must return Err, not panic"
        );
    }

    #[test]
    fn zero_bumps_returns_err() {
        let mut s = RomWorkbenchState::default();
        s.params.num_bumps = 0;
        assert!(s.run().is_err(), "num_bumps = 0 must return Err, not panic");
    }

    #[test]
    fn out_of_range_view_snapshot_is_clamped_not_panic() {
        // A view index past the last snapshot must clamp to the last column, not
        // panic.
        let mut s = RomWorkbenchState::default();
        s.params.view_snapshot = 100_000;
        let res = s.run().expect("clamped view must still run");
        assert_eq!(
            res.snapshot_index,
            s.params.num_snapshots - 1,
            "view snapshot must clamp to the last column"
        );
    }

    #[test]
    fn noise_lifts_the_spectrum_tail() {
        // With noise, the field is full-rank-ish: more significant singular
        // values than the clean low-rank field.
        let mut clean = RomWorkbenchState::default();
        clean.params.noise = 0.0;
        let mut noisy = RomWorkbenchState::default();
        noisy.params.noise = 0.2;
        let n_clean = clean.run().expect("clean").singular_values.len();
        let n_noisy = noisy.run().expect("noisy").singular_values.len();
        assert!(
            n_noisy >= n_clean,
            "noise should not reduce the significant-mode count: {n_noisy} < {n_clean}"
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
            draw_rom_workbench(app, ctx);
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
        assert!(!app.show_rom_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_rom_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rom_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rom_workbench = true;
        let res = app.rom.run().expect("run should succeed");
        app.rom.result = Some(res);
        app.rom.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_rom_workbench = true;
        // Trigger an error state (rank 0 is fail-loud in valenx-rom).
        app.rom.params.rank = 0;
        let result = app.rom.run();
        app.rom.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.rom.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_rom_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // 6 field params + 2 POD params = 8 DragValues, all exposed as
        // SpinButton nodes that MUST carry an accessible name.
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 8,
            "expected at least 8 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check specific captions are present as named accessibility nodes.
        for caption in [
            "spatial resolution n",
            "number of snapshots m",
            "number of bumps",
            "feature width \u{03C3}",
            "centre travel",
            "noise amplitude",
            "truncation rank k",
            "view snapshot index",
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
    fn numeric_controls_are_named_and_associated() {
        // Every numeric DragValue is a SpinButton and must be `labelled_by` its
        // caption (egui clears a DragValue's own Name), so an AI / screen reader
        // can find the control by its caption text. Each `labelled_by` target
        // must RESOLVE to a real named caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_rom_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 8,
            "expected the numeric controls as spin buttons, got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| {
                n.labelled_by()
                    .iter()
                    .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every DragValue's labelled_by must point at a named caption node"
        );
        for caption in ["spatial resolution n", "truncation rank k"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn degenerate_rank_shows_error_not_panic() {
        // When rank is 0 or too large the workbench must surface the error
        // in-panel, not panic.
        let mut state = RomWorkbenchState::default();
        state.params.rank = 0;
        assert!(state.run().is_err(), "rank = 0 must produce Err, not panic");
        state.params.rank = 100_000;
        assert!(
            state.run().is_err(),
            "rank > significant modes must produce Err, not panic"
        );
    }

    #[test]
    fn more_modes_lower_reconstruction_error_pin() {
        // Pinned property: more retained modes -> lower (or equal) reconstruction
        // error. (Mirrors the unit test but exercised from the UI-state struct.)
        let mut few = RomWorkbenchState::default();
        few.params.rank = 1;
        let mut many = RomWorkbenchState::default();
        many.params.rank = 7;
        let e_few = few.run().expect("few").reconstruction_error;
        let e_many = many.run().expect("many").reconstruction_error;
        assert!(
            e_many <= e_few + 1e-12,
            "more modes must not increase error: {e_many} > {e_few}"
        );
    }

    #[test]
    fn agent_bridge_rom_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "rom" }`:
        //   1. TabKind::from_id("rom") → Some(TabKind::Rom)
        //   2. set_workbench_flag(app, "rom", true) → show_rom_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup.
        assert_eq!(
            TabKind::from_id("rom"),
            Some(TabKind::Rom),
            "\"rom\" must resolve to TabKind::Rom"
        );
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("ROM"), Some(TabKind::Rom));
        assert_eq!(TabKind::from_id("  rom  "), Some(TabKind::Rom));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_rom_workbench);
        set_workbench_flag(&mut app, "rom", true);
        assert!(
            app.show_rom_workbench,
            "set_workbench_flag(\"rom\", true) must set show_rom_workbench"
        );
        set_workbench_flag(&mut app, "rom", false);
        assert!(!app.show_rom_workbench);
    }
}
