//! The right-side **UQ Workbench** panel — a native front-end over the
//! in-house `valenx-uq` crate (uncertainty quantification).
//!
//! This workbench demonstrates **forward uncertainty propagation**, **global
//! sensitivity analysis**, and **reliability analysis** on a small, fully
//! transparent model. The user describes two uncertain input random variables
//! `x1`, `x2` (each [`valenx_uq::Distribution::Normal`] or
//! [`valenx_uq::Distribution::Uniform`]), picks a model response
//! `g(x)` from a few presets (sum, product, or a linear
//! `a0 + a1·x1 + a2·x2`), sets the Monte-Carlo sample count, and clicks
//! **Run**. valenx-uq then:
//!
//! * samples the inputs with the deterministic [`valenx_uq::sampling::monte_carlo`]
//!   sampler (seeded [`valenx_uq::SplitMix64`]) and evaluates `g` over every
//!   sample into an **output sample set**;
//! * summarises that set — mean, standard deviation, and the 5/50/95 %
//!   quantiles ([`valenx_uq::statistics`]);
//! * estimates **first-order Sobol sensitivity indices** `S_i` (one per input)
//!   with the Saltelli A/B/AB design ([`valenx_uq::sensitivity::sobol_indices`]);
//! * runs **FORM** ([`valenx_uq::reliability::form`]) for the limit-state
//!   `g(x) ≤ 0`, reporting the reliability index `β` and the probability of
//!   failure `Pf = Φ(−β)`.
//!
//! Mirrors the other workbenches (`rom_workbench`, `ocean_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_uq_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"uq"` (see
//! [`crate::project_tabs::TabKind`] / [`crate::agent_commands`]).
//!
//! Two painter views are drawn: (a) a **histogram** of the output sample
//! distribution, and (b) a **bar chart** of the first-order Sobol indices (one
//! bar per input). Readouts below give the output mean / std, the failure
//! probability `Pf`, and the FORM reliability index `β`.
//!
//! Honesty: `valenx-uq` is **research / educational-grade textbook UQ** —
//! Monte-Carlo error is `O(1/√n)`, the Saltelli Sobol estimator has
//! finite-sample variance (indices may fall slightly outside `[0, 1]` at small
//! `n`), and FORM is a first-order reliability approximation (exact only for a
//! linear limit-state with Normal inputs, where `β = a0 / √(Σ aᵢ²)`). The
//! model `g` here is a deliberately simple analytic teaching example, not a
//! physical solve. Every error from `valenx-uq` surfaces verbatim — the
//! workbench never invents a number, and degenerate parameters (e.g. `0`
//! samples, a non-positive standard deviation, `lo ≥ hi`) show an in-panel
//! error, NOT a panic.

use eframe::egui;
use valenx_uq::reliability::{form, FormConfig};
use valenx_uq::sensitivity::sobol_indices;
use valenx_uq::{statistics, Distribution, FnModel, SplitMix64};

use crate::ValenxApp;

/// Fixed seed for the deterministic sampler / Sobol design so a run is
/// reproducible across machines (valenx-uq takes no `rand` dependency).
const UQ_SEED: u64 = 0x5EED_C0DE_F00D_1234;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Which family of input distribution an uncertain variable uses.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DistKind {
    /// Gaussian `N(mean, std)` — parameterised by `mean` + `std`.
    Normal,
    /// Continuous uniform `U(lo, hi)` — parameterised by `lo` + `hi`.
    Uniform,
}

/// One uncertain model input: a distribution family plus its two parameters.
///
/// `Normal` reads `p0` as the mean and `p1` as the standard deviation;
/// `Uniform` reads `p0` as the lower bound and `p1` as the upper bound. The
/// two-parameter encoding keeps the form trivial to edit and lets the same
/// drag-values serve both families (their captions change with the family).
#[derive(Clone, Copy, Debug)]
pub struct InputVar {
    /// Distribution family.
    pub kind: DistKind,
    /// First parameter — `mean` (Normal) or `lo` (Uniform).
    pub p0: f64,
    /// Second parameter — `std` (Normal) or `hi` (Uniform).
    pub p1: f64,
}

impl InputVar {
    /// Build the validated [`Distribution`] for this input, surfacing
    /// `valenx-uq`'s fail-loud constructor error verbatim (never panics).
    fn distribution(&self, label: &str) -> Result<Distribution, String> {
        match self.kind {
            DistKind::Normal => Distribution::normal(self.p0, self.p1)
                .map_err(|e| format!("input {label} (Normal): {e}")),
            DistKind::Uniform => Distribution::uniform(self.p0, self.p1)
                .map_err(|e| format!("input {label} (Uniform): {e}")),
        }
    }
}

/// The model response `g(x1, x2)` the workbench propagates uncertainty
/// through. Failure is defined as `g ≤ 0`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ModelPreset {
    /// `g = x1 + x2`.
    Sum,
    /// `g = x1 · x2`.
    Product,
    /// `g = a0 + a1·x1 + a2·x2` (the coefficients are editable in the form).
    Linear,
}

impl ModelPreset {
    /// Menu label.
    fn label(self) -> &'static str {
        match self {
            ModelPreset::Sum => "sum  g = x1 + x2",
            ModelPreset::Product => "product  g = x1 \u{00B7} x2",
            ModelPreset::Linear => "linear  g = a0 + a1\u{00B7}x1 + a2\u{00B7}x2",
        }
    }
}

/// Editable UQ inputs shown in the workbench: two random variables, the model
/// preset (+ its linear coefficients), the Monte-Carlo sample count, and the
/// failure threshold offset.
#[derive(Clone, Debug)]
pub struct UqParams {
    /// First uncertain input `x1`.
    pub x1: InputVar,
    /// Second uncertain input `x2`.
    pub x2: InputVar,
    /// The chosen response model `g`.
    pub model: ModelPreset,
    /// Linear-model intercept `a0` (only used when `model == Linear`).
    pub a0: f64,
    /// Linear-model coefficient `a1` on `x1` (only used when `model == Linear`).
    pub a1: f64,
    /// Linear-model coefficient `a2` on `x2` (only used when `model == Linear`).
    pub a2: f64,
    /// Number of Monte-Carlo samples `N` for forward propagation. Must be ≥ 2
    /// for a sample variance / Sobol estimate.
    pub n_samples: usize,
    /// Failure-threshold offset `t`: failure is declared when `g(x) − t ≤ 0`
    /// (i.e. the limit-state is `g − t`). With the default `t = 0` this is the
    /// plain `g ≤ 0` convention.
    pub threshold: f64,
}

impl Default for UqParams {
    fn default() -> Self {
        Self {
            // A modest, clearly-ordered pair so the default run is meaningful:
            // x1 has the larger variance, so for the linear model it should
            // also carry the larger first-order Sobol index.
            x1: InputVar {
                kind: DistKind::Normal,
                p0: 10.0,
                p1: 2.0,
            },
            x2: InputVar {
                kind: DistKind::Normal,
                p0: 5.0,
                p1: 1.0,
            },
            model: ModelPreset::Linear,
            a0: 1.0,
            a1: 1.0,
            a2: 1.0,
            n_samples: 4000,
            threshold: 0.0,
        }
    }
}

impl UqParams {
    /// Evaluate the chosen model `g` at an input vector `x = [x1, x2]`.
    fn eval_g(&self, x: &[f64]) -> f64 {
        match self.model {
            ModelPreset::Sum => x[0] + x[1],
            ModelPreset::Product => x[0] * x[1],
            ModelPreset::Linear => self.a0 + self.a1 * x[0] + self.a2 * x[1],
        }
    }

    /// Build a `Vec<Distribution>` for the two inputs, fail-loud.
    fn distributions(&self) -> Result<Vec<Distribution>, String> {
        Ok(vec![
            self.x1.distribution("x1")?,
            self.x2.distribution("x2")?,
        ])
    }
}

/// Parse a distribution-family name (for the agent `SetControl` bridge) into a
/// [`DistKind`]. Case-insensitive; accepts the menu words. Fail-loud on an
/// unrecognised name so a typo is a `warn` note, not a silent no-op.
fn parse_dist_kind(s: &str) -> Result<DistKind, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "normal" | "gaussian" => Ok(DistKind::Normal),
        "uniform" => Ok(DistKind::Uniform),
        other => Err(format!(
            "unknown distribution '{other}' (expected 'normal' or 'uniform')"
        )),
    }
}

/// Parse a response-model name (for the agent `SetControl` bridge) into a
/// [`ModelPreset`]. Case-insensitive; accepts the short menu words. Fail-loud.
fn parse_model_preset(s: &str) -> Result<ModelPreset, String> {
    match s.trim().to_ascii_lowercase().as_str() {
        "sum" => Ok(ModelPreset::Sum),
        "product" => Ok(ModelPreset::Product),
        "linear" => Ok(ModelPreset::Linear),
        other => Err(format!(
            "unknown model '{other}' (expected 'sum', 'product', or 'linear')"
        )),
    }
}

// ---------------------------------------------------------------------------
// Simulation result
// ---------------------------------------------------------------------------

/// Cached UQ output for the painter + readouts.
#[derive(Default, Clone)]
pub struct UqResult {
    /// The output samples `g(xᵢ)` (one per Monte-Carlo input draw), used for
    /// the histogram and the summary statistics.
    pub output_samples: Vec<f64>,
    /// Mean of the output sample set.
    pub mean: f64,
    /// Sample standard deviation of the output set (denominator `n − 1`).
    pub std: f64,
    /// 5th percentile of the output set.
    pub q05: f64,
    /// Median (50th percentile) of the output set.
    pub q50: f64,
    /// 95th percentile of the output set.
    pub q95: f64,
    /// First-order Sobol index `S_i` per input (length 2). Estimated with the
    /// Saltelli design over the **same** model.
    pub sobol_first: Vec<f64>,
    /// FORM reliability index `β = ‖u*‖` for the limit-state `g − t ≤ 0`.
    pub beta: f64,
    /// FORM probability of failure `Pf = Φ(−β)`.
    pub pf: f64,
    /// Number of HLRF iterations FORM took to converge.
    pub form_iterations: usize,
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the UQ workbench.
#[derive(Default)]
pub struct UqWorkbenchState {
    /// User-editable parameters.
    pub params: UqParams,
    /// Last successful result (populated after a successful run).
    pub result: Option<UqResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl UqWorkbenchState {
    /// The user-visible captions of every control the agent bridge can set via
    /// `SetControl` (see [`crate::agent_commands`]). Returned by
    /// `ListControls` so an agent can discover the name space. The `x1`/`x2`
    /// parameter captions are listed in their **Normal** spelling (`mean`/`std`);
    /// [`agent_set`](Self::agent_set) also accepts the **Uniform** spelling
    /// (`lo`/`hi`) for the same fields.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "x1 distribution",
            "x1 mean",
            "x1 std",
            "x2 distribution",
            "x2 mean",
            "x2 std",
            "response model g",
            "a0 (intercept)",
            "a1 (coeff on x1)",
            "a2 (coeff on x2)",
            "Monte-Carlo samples N",
            "failure threshold t",
        ]
    }

    /// Set one labelled control by its user-visible caption, for the agent
    /// `SetControl` bridge. The caption strings match exactly what the workbench
    /// form draws (and what each control is `labelled_by`), so an agent can set a
    /// parameter by the same name a user reads.
    ///
    /// Fail-loud: an unknown caption or a value of the wrong type returns
    /// `Err(String)` (the bridge turns it into a `warn` feed note) — never a
    /// panic, and no field is written on error. Numeric fields read
    /// [`AgentValue::as_f64`] / [`AgentValue::as_i64`]; the two enum captions
    /// (`x1 distribution` / `x2 distribution`, `response model g`) read
    /// [`AgentValue::as_str`] and match a small set of names. The `x1`/`x2`
    /// parameter captions accept both the Normal (`mean`/`std`) and Uniform
    /// (`lo`/`hi`) spellings since they address the same two fields (`p0`/`p1`).
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let p = &mut self.params;
        match name {
            // -- Distribution-family enums --
            "x1 distribution" => p.x1.kind = parse_dist_kind(value.as_str()?)?,
            "x2 distribution" => p.x2.kind = parse_dist_kind(value.as_str()?)?,
            // -- x1 / x2 parameters (Normal mean/std == Uniform lo/hi == p0/p1) --
            "x1 mean" | "x1 lo" => p.x1.p0 = value.as_f64()?,
            "x1 std" | "x1 hi" => p.x1.p1 = value.as_f64()?,
            "x2 mean" | "x2 lo" => p.x2.p0 = value.as_f64()?,
            "x2 std" | "x2 hi" => p.x2.p1 = value.as_f64()?,
            // -- Model enum + linear coefficients --
            "response model g" => p.model = parse_model_preset(value.as_str()?)?,
            "a0 (intercept)" => p.a0 = value.as_f64()?,
            "a1 (coeff on x1)" => p.a1 = value.as_f64()?,
            "a2 (coeff on x2)" => p.a2 = value.as_f64()?,
            // -- Sampling & limit-state --
            "Monte-Carlo samples N" => {
                let n = value.as_i64()?;
                if n < 0 {
                    return Err(format!("Monte-Carlo samples N must be >= 0, got {n}"));
                }
                p.n_samples = n as usize;
            }
            "failure threshold t" => p.threshold = value.as_f64()?,
            other => return Err(format!("unknown UQ control: {other:?}")),
        }
        Ok(())
    }

    /// The current computed-result text for the agent `ReadReadout` bridge (see
    /// [`crate::agent_commands`]). This workbench keeps its result as a structured
    /// [`UqResult`] and renders a one-line `status` summary (the `✔ mean … Pf …`
    /// line on success, or a `⚠ …` line on error) — that same `status` string is
    /// returned here. `None` when it is empty, i.e. the pipeline has not been run
    /// yet. Read-only — lets an agent read the answer back after driving a run,
    /// closing the live-driving loop.
    pub fn agent_readout(&self) -> Option<String> {
        if self.status.is_empty() {
            None
        } else {
            Some(self.status.clone())
        }
    }

    /// Run the full UQ pipeline: validate + build the input distributions,
    /// Monte-Carlo sample them, evaluate `g` to an output sample set, summarise
    /// it (mean / std / quantiles), estimate first-order Sobol indices
    /// (Saltelli), and run FORM for the limit-state `g − t ≤ 0`.
    ///
    /// Every failure is returned as an `Err(String)` — no panics, no invented
    /// numbers. Degenerate inputs (`n_samples < 2`, a non-positive Normal
    /// `std`, a `Uniform` with `lo ≥ hi`, a constant response that carries no
    /// apportionable Sobol variance, or a FORM that does not converge) surface
    /// `valenx-uq`'s own error / `None` verbatim.
    pub fn run(&self) -> Result<UqResult, String> {
        let p = &self.params;

        if p.n_samples < 2 {
            return Err(format!(
                "Monte-Carlo needs >= 2 samples for a variance / Sobol estimate, got {}",
                p.n_samples
            ));
        }

        // --- Input distributions (fail-loud on bad params) ------------------
        let dists = p.distributions()?;

        // --- Forward propagation: sample, evaluate g, summarise ------------
        let mut rng = SplitMix64::new(UQ_SEED);
        let inputs = valenx_uq::sampling::monte_carlo(p.n_samples, &dists, &mut rng);
        let output_samples: Vec<f64> = inputs.iter().map(|x| p.eval_g(x)).collect();

        let mean = statistics::mean(&output_samples)
            .ok_or_else(|| "output sample set is empty (cannot take a mean)".to_string())?;
        let std = statistics::std(&output_samples)
            .ok_or_else(|| "need >= 2 output samples for a standard deviation".to_string())?;
        let q05 = statistics::percentile(&output_samples, 5.0).map_err(|e| e.to_string())?;
        let q50 = statistics::median(&output_samples).map_err(|e| e.to_string())?;
        let q95 = statistics::percentile(&output_samples, 95.0).map_err(|e| e.to_string())?;

        // --- First-order Sobol indices (Saltelli A/B/AB) -------------------
        // Wrap g behind the Model trait. The Saltelli base sample size is
        // bounded so the (d+2)·N evaluation budget stays interactive even at a
        // large forward-propagation N.
        let g_closure = {
            let p = p.clone();
            move |x: &[f64]| vec![p.eval_g(x)]
        };
        let model = FnModel::new(2, 1, g_closure);
        let n_base = p.n_samples.clamp(2, 8192);
        let mut sobol_rng = SplitMix64::new(UQ_SEED ^ 0xA5A5_A5A5_A5A5_A5A5);
        let sobol = sobol_indices(&model, &dists, n_base, 0, &mut sobol_rng).ok_or_else(|| {
            "Sobol analysis is ill-posed (constant response → zero output variance, \
             or too few samples). Try a non-degenerate model / more samples."
                .to_string()
        })?;
        let sobol_first = sobol.first_order;

        // --- FORM reliability for the limit-state g - t <= 0 ---------------
        let threshold = p.threshold;
        let limit_state = {
            let p = p.clone();
            move |x: &[f64]| p.eval_g(x) - threshold
        };
        let form_r =
            form(limit_state, &dists, &FormConfig::default()).map_err(|e| e.to_string())?;

        Ok(UqResult {
            output_samples,
            mean,
            std,
            q05,
            q50,
            q95,
            sobol_first,
            beta: form_r.beta,
            pf: form_r.pf,
            form_iterations: form_r.iterations,
        })
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the UQ workbench. A no-op unless toggled on via View → UQ.
///
/// Mirrors [`crate::rom_workbench::draw_rom_workbench`].
pub fn draw_uq_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_uq_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_uq_workbench",
        "Uncertainty quantification",
        uq_workbench_body,
    );
    if close {
        app.show_uq_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn uq_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Uncertainty quantification \u{2014} forward Monte-Carlo propagation, first-order \
             Sobol sensitivity (Saltelli), and FORM reliability (Pf) for the limit-state \
             g \u{2264} 0 \u{00B7} valenx-uq  [research / educational \u{2014} MC error O(1/\u{221A}n); \
             FORM is first-order]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.uq;
        let p = &mut s.params;

        ui.label(egui::RichText::new("Input random variables").strong());
        egui::Grid::new("uq_input_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                input_var_rows(ui, "x1", &mut p.x1);
                input_var_rows(ui, "x2", &mut p.x2);
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Model g(x)").strong());
        egui::Grid::new("uq_model_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("response model g");
                egui::ComboBox::from_id_source("uq_model_combo")
                    .selected_text(p.model.label())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut p.model,
                            ModelPreset::Sum,
                            ModelPreset::Sum.label(),
                        );
                        ui.selectable_value(
                            &mut p.model,
                            ModelPreset::Product,
                            ModelPreset::Product.label(),
                        );
                        ui.selectable_value(
                            &mut p.model,
                            ModelPreset::Linear,
                            ModelPreset::Linear.label(),
                        );
                    })
                    .response
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "The model response g(x1, x2). Failure (for FORM / Pf) is g <= threshold.",
                    );
                ui.end_row();

                // Linear coefficients only matter for the linear preset; show
                // them always (greyed when unused) so the form layout is stable
                // and the controls keep stable accessible names.
                let linear = p.model == ModelPreset::Linear;
                ui.add_enabled_ui(linear, |ui| {
                    let lbl = ui.label("a0 (intercept)");
                    ui.add(egui::DragValue::new(&mut p.a0).speed(0.1))
                        .labelled_by(lbl.id)
                        .on_hover_text("Linear-model intercept a0 in g = a0 + a1*x1 + a2*x2.");
                });
                ui.end_row();
                ui.add_enabled_ui(linear, |ui| {
                    let lbl = ui.label("a1 (coeff on x1)");
                    ui.add(egui::DragValue::new(&mut p.a1).speed(0.1))
                        .labelled_by(lbl.id)
                        .on_hover_text("Linear-model coefficient a1 on x1.");
                });
                ui.end_row();
                ui.add_enabled_ui(linear, |ui| {
                    let lbl = ui.label("a2 (coeff on x2)");
                    ui.add(egui::DragValue::new(&mut p.a2).speed(0.1))
                        .labelled_by(lbl.id)
                        .on_hover_text("Linear-model coefficient a2 on x2.");
                });
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Sampling & limit-state").strong());
        egui::Grid::new("uq_run_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("Monte-Carlo samples N");
                ui.add(
                    egui::DragValue::new(&mut p.n_samples)
                        .speed(50)
                        .range(0..=2_000_000),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Number of input draws propagated through g. Must be >= 2; \
                     Monte-Carlo error shrinks as O(1/sqrt(N)).",
                );
                ui.end_row();

                let lbl = ui.label("failure threshold t");
                ui.add(egui::DragValue::new(&mut p.threshold).speed(0.1))
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Limit-state offset: failure is declared when g(x) - t <= 0. \
                         Leave at 0 for the plain g <= 0 convention.",
                    );
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Sample the inputs, evaluate g, and compute statistics, Sobol indices, \
                     and FORM reliability.",
                )
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
    let s = &app.uq;
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
    draw_uq_viz(s, ui);
}

/// One input variable's three rows in the parameter grid: a distribution-family
/// combo and the two parameter drag-values whose captions follow the family.
fn input_var_rows(ui: &mut egui::Ui, name: &str, v: &mut InputVar) {
    let lbl = ui.label(format!("{name} distribution"));
    egui::ComboBox::from_id_source(format!("uq_{name}_dist"))
        .selected_text(match v.kind {
            DistKind::Normal => "Normal (mean, std)",
            DistKind::Uniform => "Uniform (lo, hi)",
        })
        .show_ui(ui, |ui| {
            ui.selectable_value(&mut v.kind, DistKind::Normal, "Normal (mean, std)");
            ui.selectable_value(&mut v.kind, DistKind::Uniform, "Uniform (lo, hi)");
        })
        .response
        .labelled_by(lbl.id)
        .on_hover_text("Distribution family for this uncertain input.");
    ui.end_row();

    let (p0_cap, p1_cap, p0_hint, p1_hint) = match v.kind {
        DistKind::Normal => (
            format!("{name} mean"),
            format!("{name} std"),
            "Mean of the Normal input.",
            "Standard deviation of the Normal input. Must be > 0.",
        ),
        DistKind::Uniform => (
            format!("{name} lo"),
            format!("{name} hi"),
            "Lower bound of the Uniform input.",
            "Upper bound of the Uniform input. Must be > lo.",
        ),
    };

    let lbl = ui.label(p0_cap);
    ui.add(egui::DragValue::new(&mut v.p0).speed(0.1))
        .labelled_by(lbl.id)
        .on_hover_text(p0_hint);
    ui.end_row();

    let lbl = ui.label(p1_cap);
    ui.add(egui::DragValue::new(&mut v.p1).speed(0.1))
        .labelled_by(lbl.id)
        .on_hover_text(p1_hint);
    ui.end_row();
}

/// Run the pipeline and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.uq;
    match s.run() {
        Ok(res) => {
            s.status = format!(
                "\u{2714} mean {:.4} \u{00B7} std {:.4} \u{00B7} Pf {:.3e} \u{00B7} \u{03B2} {:.4} \
                 \u{00B7} S1 {:.3} / S2 {:.3}",
                res.mean,
                res.std,
                res.pf,
                res.beta,
                res.sobol_first.first().copied().unwrap_or(f64::NAN),
                res.sobol_first.get(1).copied().unwrap_or(f64::NAN),
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
// 2-D visualisation (output histogram + Sobol-index bar chart)
// ---------------------------------------------------------------------------

fn draw_uq_viz(s: &UqWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to propagate uncertainty and visualise the output distribution",
            )
            .weak(),
        );
        return;
    };

    draw_output_histogram(res, ui);
    ui.add_space(8.0);
    draw_sobol_bars(res, ui);

    // Readouts grid below the painters.
    ui.add_space(6.0);
    egui::Grid::new("uq_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(ui, "output mean", format!("{:.5}", res.mean));
            row(ui, "output std", format!("{:.5}", res.std));
            row(
                ui,
                "output 5 / 50 / 95 %",
                format!("{:.4} / {:.4} / {:.4}", res.q05, res.q50, res.q95),
            );
            row(
                ui,
                "Sobol S1 (x1) / S2 (x2)",
                format!(
                    "{:.4} / {:.4}",
                    res.sobol_first.first().copied().unwrap_or(f64::NAN),
                    res.sobol_first.get(1).copied().unwrap_or(f64::NAN),
                ),
            );
            row(
                ui,
                "FORM reliability index \u{03B2}",
                format!("{:.5}", res.beta),
            );
            row(ui, "probability of failure Pf", format!("{:.5e}", res.pf));
            row(
                ui,
                "FORM HLRF iterations",
                format!("{}", res.form_iterations),
            );
            row(
                ui,
                "output samples N",
                format!("{}", res.output_samples.len()),
            );
        });
}

/// View (a): a histogram of the output sample distribution.
fn draw_output_histogram(res: &UqResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("Output distribution histogram").strong());
    ui.label(
        egui::RichText::new(
            "teal bars = g(x) sample density \u{00B7} green dashed = mean \u{00B7} \
             grey dashed = failure threshold (g = 0)",
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

    if res.output_samples.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "too few samples to histogram",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    // Sample range. A degenerate (constant) sample set would give a zero-width
    // range; pad it so the single bar is still drawn.
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in &res.output_samples {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    if !(lo.is_finite() && hi.is_finite()) {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "non-finite samples",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }
    if (hi - lo).abs() < 1e-12 {
        lo -= 0.5;
        hi += 0.5;
    }
    let span = hi - lo;

    let margin = 16.0_f32;
    let inner = rect.shrink(margin);

    // Bin the samples (Sturges-ish bin count, bounded).
    let n_bins = ((res.output_samples.len() as f64).sqrt().ceil() as usize).clamp(8, 48);
    let mut bins = vec![0_usize; n_bins];
    for &v in &res.output_samples {
        if !v.is_finite() {
            continue;
        }
        let frac = ((v - lo) / span).clamp(0.0, 1.0);
        let b = ((frac * n_bins as f64) as usize).min(n_bins - 1);
        bins[b] += 1;
    }
    let max_count = bins.iter().copied().max().unwrap_or(1).max(1);

    let bar_w = inner.width() / n_bins as f32;
    for (i, &c) in bins.iter().enumerate() {
        let h = (c as f32 / max_count as f32) * inner.height();
        let x0 = inner.left() + i as f32 * bar_w;
        let bar = egui::Rect::from_min_max(
            egui::pos2(x0 + 0.5, inner.bottom() - h),
            egui::pos2(x0 + bar_w - 0.5, inner.bottom()),
        );
        painter.rect_filled(bar, 0.0, egui::Color32::from_rgb(70, 170, 180));
    }

    // Map an output value to an x pixel inside the inner rect.
    let val_to_x = |v: f64| -> f32 {
        let frac = ((v - lo) / span).clamp(0.0, 1.0) as f32;
        inner.left() + frac * inner.width()
    };

    // Dashed vertical line at the mean.
    dashed_vline(
        &painter,
        val_to_x(res.mean),
        inner.top(),
        inner.bottom(),
        egui::Color32::from_rgb(110, 200, 130),
    );

    // Dashed vertical line at the failure threshold g = 0, if it's in range.
    if 0.0 >= lo && 0.0 <= hi {
        dashed_vline(
            &painter,
            val_to_x(0.0),
            inner.top(),
            inner.bottom(),
            egui::Color32::from_gray(150),
        );
    }
}

/// View (b): a bar chart of the first-order Sobol indices, one bar per input.
fn draw_sobol_bars(res: &UqResult, ui: &mut egui::Ui) {
    ui.label(egui::RichText::new("First-order Sobol indices S_i").strong());
    ui.label(
        egui::RichText::new(
            "amber bars = fraction of output variance from each input alone \
             (Saltelli estimate; may be slightly outside [0,1] at small N)",
        )
        .weak()
        .small(),
    );

    let available = ui.available_size();
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(available.x.min(460.0), 160.0),
        egui::Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.sobol_first.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no Sobol indices",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let margin = 22.0_f32;
    let inner = rect.shrink(margin);

    // Map [0, 1] of the index to the full bar height; clamp the (estimator)
    // overshoot/undershoot so a slightly-out-of-range index still draws sanely.
    let count = res.sobol_first.len();
    let slot_w = inner.width() / count as f32;
    let labels = ["x1", "x2", "x3"];

    // Zero baseline.
    painter.line_segment(
        [
            egui::pos2(inner.left(), inner.bottom()),
            egui::pos2(inner.right(), inner.bottom()),
        ],
        egui::Stroke::new(1.0, egui::Color32::from_gray(70)),
    );

    for (i, &s) in res.sobol_first.iter().enumerate() {
        let frac = (s.clamp(0.0, 1.0)) as f32;
        let h = frac * inner.height();
        let cx = inner.left() + (i as f32 + 0.5) * slot_w;
        let half = (slot_w * 0.30).max(4.0);
        let bar = egui::Rect::from_min_max(
            egui::pos2(cx - half, inner.bottom() - h),
            egui::pos2(cx + half, inner.bottom()),
        );
        painter.rect_filled(bar, 0.0, egui::Color32::from_rgb(230, 180, 70));

        // Input label under the bar.
        let name = labels.get(i).copied().unwrap_or("x");
        painter.text(
            egui::pos2(cx, inner.bottom() + 4.0),
            egui::Align2::CENTER_TOP,
            name,
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(180),
        );
        // Numeric value above the bar.
        painter.text(
            egui::pos2(cx, inner.bottom() - h - 2.0),
            egui::Align2::CENTER_BOTTOM,
            format!("{s:.3}"),
            egui::FontId::monospace(11.0),
            egui::Color32::from_rgb(240, 220, 170),
        );
    }
}

/// Draw a short-dashed vertical line between `y0` and `y1` at pixel `x`.
fn dashed_vline(painter: &egui::Painter, x: f32, y0: f32, y1: f32, color: egui::Color32) {
    let mut y = y0;
    while y < y1 {
        let y2 = (y + 5.0).min(y1);
        painter.line_segment(
            [egui::pos2(x, y), egui::pos2(x, y2)],
            egui::Stroke::new(1.0, color),
        );
        y += 9.0;
    }
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring rom_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_run_succeeds_and_is_populated() {
        let s = UqWorkbenchState::default();
        let res = s.run().expect("default UQ run should succeed");
        assert_eq!(
            res.output_samples.len(),
            s.params.n_samples,
            "one output sample per input draw"
        );
        assert!(res.mean.is_finite(), "mean must be finite");
        assert!(
            res.std.is_finite() && res.std >= 0.0,
            "std must be finite >= 0"
        );
        assert_eq!(res.sobol_first.len(), 2, "two inputs -> two Sobol indices");
        assert!(res.beta.is_finite() && res.beta >= 0.0, "beta finite >= 0");
        assert!(
            (0.0..=1.0).contains(&res.pf),
            "Pf must be in [0,1], got {}",
            res.pf
        );
        // Quantiles are ordered.
        assert!(
            res.q05 <= res.q50 + 1e-9 && res.q50 <= res.q95 + 1e-9,
            "quantiles must be ordered: {} {} {}",
            res.q05,
            res.q50,
            res.q95
        );
    }

    #[test]
    fn linear_mean_matches_closed_form() {
        // For g = a0 + a1*x1 + a2*x2 with E[x1]=m1, E[x2]=m2, the output mean
        // is a0 + a1*m1 + a2*m2 (up to Monte-Carlo error).
        let s = UqWorkbenchState::default();
        let p = &s.params;
        let expected = p.a0 + p.a1 * p.x1.p0 + p.a2 * p.x2.p0;
        let res = s.run().expect("run");
        assert!(
            (res.mean - expected).abs() < 0.2,
            "output mean {} should be near closed-form {expected}",
            res.mean
        );
    }

    #[test]
    fn form_beta_for_linear_limit_state_matches_a0_over_norm() {
        // Pinned property: for a LINEAR limit-state g = a0 + a1*x1 + a2*x2 with
        // standard-normal inputs and threshold 0, FORM is exact:
        //   beta = a0 / sqrt(a1^2 + a2^2).
        let mut s = UqWorkbenchState::default();
        s.params.x1 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 1.0,
        };
        s.params.x2 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 1.0,
        };
        s.params.model = ModelPreset::Linear;
        s.params.a0 = 3.0;
        s.params.a1 = 2.0;
        s.params.a2 = 1.0;
        s.params.threshold = 0.0;

        let res = s.run().expect("run");
        let beta_exact = 3.0_f64 / (2.0_f64 * 2.0 + 1.0).sqrt();
        assert!(
            (res.beta - beta_exact).abs() < 1e-5,
            "FORM beta {} should equal a0/sqrt(a1^2+a2^2) = {beta_exact}",
            res.beta
        );
    }

    #[test]
    fn larger_variance_input_has_larger_sobol_index() {
        // Pinned property: for g = a1*x1 + a2*x2 (linear, independent inputs),
        // the first-order Sobol index of an input is proportional to its
        // variance contribution (a_i^2 * Var(x_i)). With equal coefficients,
        // the input with the larger variance must carry the larger S_i.
        let mut s = UqWorkbenchState::default();
        s.params.model = ModelPreset::Linear;
        s.params.a0 = 0.0;
        s.params.a1 = 1.0;
        s.params.a2 = 1.0;
        // x1 std 3 (var 9) >> x2 std 1 (var 1).
        s.params.x1 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 3.0,
        };
        s.params.x2 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 1.0,
        };
        s.params.n_samples = 4096;

        let res = s.run().expect("run");
        let s1 = res.sobol_first[0];
        let s2 = res.sobol_first[1];
        assert!(
            s1 > s2,
            "x1 (var 9) must have a larger Sobol index than x2 (var 1): S1={s1}, S2={s2}"
        );
    }

    // ---- degenerate-param tests — must return Err, NOT panic ----

    #[test]
    fn zero_samples_returns_err() {
        let mut s = UqWorkbenchState::default();
        s.params.n_samples = 0;
        assert!(s.run().is_err(), "n_samples = 0 must return Err, not panic");
    }

    #[test]
    fn nonpositive_std_returns_err() {
        let mut s = UqWorkbenchState::default();
        s.params.x1.kind = DistKind::Normal;
        s.params.x1.p1 = 0.0; // std <= 0
        assert!(s.run().is_err(), "std <= 0 must return Err, not panic");

        s.params.x1.p1 = -1.0;
        assert!(s.run().is_err(), "negative std must return Err, not panic");
    }

    #[test]
    fn uniform_lo_ge_hi_returns_err() {
        let mut s = UqWorkbenchState::default();
        s.params.x2.kind = DistKind::Uniform;
        s.params.x2.p0 = 5.0; // lo
        s.params.x2.p1 = 5.0; // hi == lo
        assert!(
            s.run().is_err(),
            "uniform lo >= hi must return Err, not panic"
        );
    }

    #[test]
    fn uniform_inputs_run_without_error() {
        // Sanity: a valid Uniform/Uniform pair propagates fine.
        let mut s = UqWorkbenchState::default();
        s.params.x1 = InputVar {
            kind: DistKind::Uniform,
            p0: 0.0,
            p1: 10.0,
        };
        s.params.x2 = InputVar {
            kind: DistKind::Uniform,
            p0: 1.0,
            p1: 3.0,
        };
        s.params.model = ModelPreset::Sum;
        let res = s.run().expect("uniform run should succeed");
        assert_eq!(res.sobol_first.len(), 2);
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
            draw_uq_workbench(app, ctx);
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
        assert!(!app.show_uq_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_uq_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_uq_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_uq_workbench = true;
        let res = app.uq.run().expect("run should succeed");
        app.uq.result = Some(res);
        app.uq.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_uq_workbench = true;
        // Trigger an error state (0 samples is fail-loud in run()).
        app.uq.params.n_samples = 0;
        let result = app.uq.run();
        app.uq.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.uq.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_uq_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // Numeric DragValues: 4 input params (x1.p0/p1, x2.p0/p1) + 3 linear
        // coeffs (a0/a1/a2) + 2 run params (N, threshold) = 9 SpinButtons, all
        // of which MUST carry an accessible name (be labelled_by a caption).
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 9,
            "expected at least 9 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check specific captions are present as named accessibility nodes
        // (default x1/x2 are Normal, so the captions are mean/std).
        for caption in [
            "x1 distribution",
            "x1 mean",
            "x1 std",
            "x2 distribution",
            "x2 mean",
            "x2 std",
            "response model g",
            "a0 (intercept)",
            "a1 (coeff on x1)",
            "a2 (coeff on x2)",
            "Monte-Carlo samples N",
            "failure threshold t",
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
        // can find the control by caption text. Each `labelled_by` target must
        // RESOLVE to a real named caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_uq_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 9,
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
        for caption in ["x1 mean", "Monte-Carlo samples N"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn degenerate_params_show_error_not_panic() {
        // When the sample count is 0 or a std is non-positive the workbench must
        // surface the error in-panel, not panic.
        let mut state = UqWorkbenchState::default();
        state.params.n_samples = 0;
        assert!(
            state.run().is_err(),
            "n_samples = 0 must produce Err, not panic"
        );
        state.params.n_samples = 4000;
        state.params.x1.p1 = 0.0; // std <= 0
        assert!(state.run().is_err(), "std <= 0 must produce Err, not panic");
    }

    #[test]
    fn form_beta_linear_pin_from_ui_state() {
        // Mirror of the unit pin, exercised from the UI-state struct: FORM beta
        // for a linear limit-state with standard-normal inputs equals
        // a0/sqrt(a1^2 + a2^2).
        let mut s = UqWorkbenchState::default();
        s.params.x1 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 1.0,
        };
        s.params.x2 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 1.0,
        };
        s.params.model = ModelPreset::Linear;
        s.params.a0 = 4.0;
        s.params.a1 = 3.0;
        s.params.a2 = 0.0;
        let res = s.run().expect("run");
        let beta_exact = 4.0_f64 / (3.0_f64 * 3.0).sqrt(); // = 4/3
        assert!(
            (res.beta - beta_exact).abs() < 1e-5,
            "beta {} should equal a0/sqrt(a1^2+a2^2) = {beta_exact}",
            res.beta
        );
    }

    #[test]
    fn larger_variance_input_has_larger_sobol_index_pin() {
        // Pinned property (from UI state): for g = a1*x1 + a2*x2 with equal
        // coefficients, the input with the larger variance has the larger
        // first-order Sobol index.
        let mut s = UqWorkbenchState::default();
        s.params.model = ModelPreset::Linear;
        s.params.a0 = 0.0;
        s.params.a1 = 1.0;
        s.params.a2 = 1.0;
        s.params.x1 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 4.0, // var 16
        };
        s.params.x2 = InputVar {
            kind: DistKind::Normal,
            p0: 0.0,
            p1: 1.0, // var 1
        };
        s.params.n_samples = 4096;
        let res = s.run().expect("run");
        assert!(
            res.sobol_first[0] > res.sobol_first[1],
            "x1 (var 16) Sobol index {} must exceed x2 (var 1) {}",
            res.sobol_first[0],
            res.sobol_first[1]
        );
    }

    #[test]
    fn agent_bridge_uq_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "uq" }`:
        //   1. TabKind::from_id("uq") -> Some(TabKind::Uq)
        //   2. set_workbench_flag(app, "uq", true) -> show_uq_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup.
        assert_eq!(
            TabKind::from_id("uq"),
            Some(TabKind::Uq),
            "\"uq\" must resolve to TabKind::Uq"
        );
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("UQ"), Some(TabKind::Uq));
        assert_eq!(TabKind::from_id("  uq  "), Some(TabKind::Uq));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_uq_workbench);
        set_workbench_flag(&mut app, "uq", true);
        assert!(
            app.show_uq_workbench,
            "set_workbench_flag(\"uq\", true) must set show_uq_workbench"
        );
        set_workbench_flag(&mut app, "uq", false);
        assert!(!app.show_uq_workbench);
    }
}
