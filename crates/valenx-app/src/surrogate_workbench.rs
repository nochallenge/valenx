//! The right-side **Surrogate Model** workbench — a native, in-house **machine
//! -learning surrogate** (a "response-surface emulator") for an expensive solver.
//!
//! ## The idea
//!
//! Many valenx solvers (FEM, CFD, …) are slow enough that you cannot re-run them
//! while *dragging a slider*. A **surrogate model** fixes that: you sample the
//! true solver at `N` design points, **train a small neural net** on those
//! `(inputs -> output)` pairs once, and from then on the net predicts the output
//! **instantly (microseconds)** for any new input. That is the
//! slider-moves-and-stress-updates-instantly experience — the surrogate stands in
//! for the solver inside the interactive loop, and you re-validate against the
//! true solver only when it matters.
//!
//! ## What V1 ships
//!
//! The ground-truth function is a **cantilever-beam tip deflection**, which has a
//! known closed form so the surrogate's accuracy is *validatable* against truth:
//!
//! ```text
//! delta(x0, x1) = P * L^3 / (3 * E * I)
//! ```
//!
//! parameterised by two normalised inputs on `[0, 1]` — a **load** factor and a
//! **length** factor — mapped to physical `P` and `L` over fixed ranges (with `E`
//! and `I` held constant), so it is a genuine, smooth, nonlinear 2-input response
//! surface (cubic in length, linear in load). See [`Truth::deflection`].
//!
//! The surrogate is a tiny **in-house multilayer perceptron** (`Mlp`):
//! `2 -> H -> H -> 1` dense layers with **ReLU** hidden activations and a linear
//! output, trained by **full-batch gradient descent with Adam** on **mean-squared
//! error**. Inputs are fed in normalised `[0,1]` space and the target is
//! standardised (zero mean / unit variance) so the net trains stably; the
//! standardisation is inverted on predict so the reported value is in physical
//! units (metres). No external DL framework is pulled in — the net is a few
//! hundred lines of plain `f64` arithmetic, deterministic under a seeded RNG, and
//! fully unit-testable headless. (`burn`/`candle` were considered per crate-first
//! but rejected: their default backends risk a `wgpu` version clash with eframe's
//! pinned wgpu/egui-0.28, and a 2->H->H->1 net does not need a tensor framework.)
//!
//! Training reports **train and test MSE** on a held-out split; a unit test
//! asserts the test MSE (in normalised target units) is small — i.e. the net
//! actually learned the function.
//!
//! ## AI-drivable surface (the #1 standing release gate)
//!
//! Mirrors every other workbench: a [`crate::workbench_chrome::workbench_shell`]
//! panel gated on [`crate::ValenxApp::show_surrogate_workbench`], toggled from the
//! View menu and openable by the agent bridge under the id `"surrogate"`. The
//! bridge can:
//! - set the two prediction-input sliders, the training-sample count, the hidden
//!   width, the epoch count and the test-split fraction through labelled controls
//!   (`agent_set` / `agent_control_names`);
//! - read a status line (`agent_readout`): the train / test MSE and the current
//!   surrogate-vs-true prediction at the input sliders;
//! - fire training via the `RunCommand` id `surrogate.train`, which runs the SAME
//!   sample-then-train path the in-panel **Train** button calls.

use eframe::egui;
use egui_plot::{Legend, Line, Plot, PlotPoints};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Deterministic RNG (tiny SplitMix64 — no `rand` dependency needed here)
// ---------------------------------------------------------------------------

/// A minimal deterministic PRNG (SplitMix64). Used for the train/test sampling
/// and the MLP weight init so a given seed yields a reproducible run (the tests
/// rely on this; the UI reseeds per Train so results are stable run-to-run).
#[derive(Clone, Debug)]
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
        }
    }

    /// Next `u64` (SplitMix64).
    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform `f64` in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        // Top 53 bits -> [0,1).
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Standard-normal `f64` via Box-Muller (for weight init).
    fn normal(&mut self) -> f64 {
        // Guard u1 away from 0 so ln() is finite.
        let u1 = self.unit().max(1e-12);
        let u2 = self.unit();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }
}

// ---------------------------------------------------------------------------
// Ground-truth function (the "expensive solver" stand-in, known closed form)
// ---------------------------------------------------------------------------

/// The ground-truth response surface the surrogate learns: **cantilever-beam tip
/// deflection** `delta = P*L^3/(3*E*I)`, as a function of two normalised inputs
/// `x0, x1` in `[0,1]` (a *load* factor and a *length* factor). `E` (Young's
/// modulus) and `I` (second moment of area) are fixed; `P` and `L` sweep over
/// fixed physical ranges. Smooth, nonlinear (cubic in `L`, linear in `P`), and —
/// crucially — *closed-form*, so the surrogate's error against truth is exactly
/// measurable.
#[derive(Clone, Copy, Debug)]
pub struct Truth;

impl Truth {
    /// Fixed Young's modulus `E` (Pa) — structural steel.
    const E: f64 = 200e9;
    /// Fixed second moment of area `I` (m^4) — a slender rectangular section.
    const I: f64 = 8.0e-6;
    /// Load `P` (N) range mapped from `x0` in `[0,1]`: `100 .. 1100 N`.
    const P_MIN: f64 = 100.0;
    const P_SPAN: f64 = 1000.0;
    /// Length `L` (m) range mapped from `x1` in `[0,1]`: `0.5 .. 2.5 m`.
    const L_MIN: f64 = 0.5;
    const L_SPAN: f64 = 2.0;

    /// Map the normalised load factor `x0` to physical load `P` (N).
    pub fn load(x0: f64) -> f64 {
        Self::P_MIN + Self::P_SPAN * x0.clamp(0.0, 1.0)
    }

    /// Map the normalised length factor `x1` to physical length `L` (m).
    pub fn length(x1: f64) -> f64 {
        Self::L_MIN + Self::L_SPAN * x1.clamp(0.0, 1.0)
    }

    /// The true tip deflection `delta` (m) for normalised inputs `x0` (load),
    /// `x1` (length). This is the "solver" the surrogate emulates.
    pub fn deflection(x0: f64, x1: f64) -> f64 {
        let p = Self::load(x0);
        let l = Self::length(x1);
        p * l * l * l / (3.0 * Self::E * Self::I)
    }
}

// ---------------------------------------------------------------------------
// Tiny in-house MLP: 2 -> H -> H -> 1, ReLU hidden, linear out, Adam + MSE
// ---------------------------------------------------------------------------

/// One fully-connected layer `y = W*x + b` with optional ReLU, plus the Adam
/// moment buffers for its parameters. Row-major weights: `w[o][i]` is the weight
/// from input `i` to output `o`.
#[derive(Clone, Debug)]
struct Layer {
    /// `out x in` weights.
    w: Vec<Vec<f64>>,
    /// `out` biases.
    b: Vec<f64>,
    /// Apply ReLU to this layer's output (hidden layers `true`, output `false`).
    relu: bool,
    // --- Adam state (same shapes as w / b) ---
    mw: Vec<Vec<f64>>,
    vw: Vec<Vec<f64>>,
    mb: Vec<f64>,
    vb: Vec<f64>,
}

impl Layer {
    /// He-initialise an `in -> out` layer (good for ReLU); biases zero.
    fn new(n_in: usize, n_out: usize, relu: bool, rng: &mut Rng) -> Self {
        let scale = (2.0 / n_in as f64).sqrt();
        let w: Vec<Vec<f64>> = (0..n_out)
            .map(|_| (0..n_in).map(|_| rng.normal() * scale).collect())
            .collect();
        Self {
            mw: vec![vec![0.0; n_in]; n_out],
            vw: vec![vec![0.0; n_in]; n_out],
            mb: vec![0.0; n_out],
            vb: vec![0.0; n_out],
            w,
            b: vec![0.0; n_out],
            relu,
        }
    }

    /// Forward pass: returns the **pre-activation** `z` and the **activation**
    /// `a` (with ReLU applied iff `self.relu`). Both are needed for backprop.
    fn forward(&self, x: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let mut z = self.b.clone();
        for (o, zo) in z.iter_mut().enumerate() {
            let wo = &self.w[o];
            let mut acc = *zo;
            for (i, &xi) in x.iter().enumerate() {
                acc += wo[i] * xi;
            }
            *zo = acc;
        }
        let a = if self.relu {
            z.iter().map(|v| v.max(0.0)).collect()
        } else {
            z.clone()
        };
        (z, a)
    }
}

/// The surrogate network: a stack of [`Layer`]s and the Adam timestep counter.
/// Predicts in **standardised target space**; the owning [`SurrogateModel`]
/// inverts the standardisation to physical units.
#[derive(Clone, Debug)]
struct Mlp {
    layers: Vec<Layer>,
    /// Adam timestep (for bias-correction).
    t: u64,
}

impl Mlp {
    /// Build a `2 -> hidden -> hidden -> 1` network with the given hidden width.
    fn new(hidden: usize, rng: &mut Rng) -> Self {
        let hidden = hidden.max(1);
        Self {
            layers: vec![
                Layer::new(2, hidden, true, rng),
                Layer::new(hidden, hidden, true, rng),
                Layer::new(hidden, 1, false, rng),
            ],
            t: 0,
        }
    }

    /// Predict the (standardised) scalar output for a 2-vector input.
    fn predict(&self, x: &[f64]) -> f64 {
        let mut a = x.to_vec();
        for l in &self.layers {
            a = l.forward(&a).1;
        }
        a[0]
    }

    /// One full-batch Adam step on MSE over `(x, y)` standardised pairs; returns
    /// the batch MSE *before* the step. Vanilla backprop through the dense + ReLU
    /// stack — kept explicit (no autodiff dep) since the graph is fixed.
    fn train_step(&mut self, xs: &[[f64; 2]], ys: &[f64], lr: f64) -> f64 {
        let n = xs.len().max(1);
        let nl = self.layers.len();

        // Accumulated parameter gradients (same shapes as each layer's w / b).
        let mut gw: Vec<Vec<Vec<f64>>> = self
            .layers
            .iter()
            .map(|l| vec![vec![0.0; l.w[0].len()]; l.w.len()])
            .collect();
        let mut gb: Vec<Vec<f64>> = self.layers.iter().map(|l| vec![0.0; l.b.len()]).collect();

        let mut mse = 0.0;
        for (sample, x) in xs.iter().enumerate() {
            // ---- forward, caching per-layer (input, pre-activation) ----
            let mut acts: Vec<Vec<f64>> = Vec::with_capacity(nl + 1);
            let mut zs: Vec<Vec<f64>> = Vec::with_capacity(nl);
            acts.push(x.to_vec());
            for l in &self.layers {
                let (z, a) = l.forward(acts.last().unwrap());
                zs.push(z);
                acts.push(a);
            }
            let pred = acts[nl][0];
            let err = pred - ys[sample];
            mse += err * err;

            // ---- backward ----
            // dL/d(output activation) for MSE = 2(pred - y); fold the 1/n here.
            let mut delta = vec![2.0 * err / n as f64];
            for li in (0..nl).rev() {
                let layer = &self.layers[li];
                let z = &zs[li];
                // Apply activation derivative for hidden ReLU layers.
                if layer.relu {
                    for (d, &zv) in delta.iter_mut().zip(z.iter()) {
                        if zv <= 0.0 {
                            *d = 0.0;
                        }
                    }
                }
                let a_in = &acts[li];
                // Param grads for this layer (`delta` has one entry per output o).
                for (o, &d) in delta.iter().enumerate() {
                    gb[li][o] += d;
                    let go = &mut gw[li][o];
                    for (i, &ai) in a_in.iter().enumerate() {
                        go[i] += d * ai;
                    }
                }
                // Propagate to the previous layer's activation (if any).
                if li > 0 {
                    let mut prev = vec![0.0; a_in.len()];
                    for (o, &d) in delta.iter().enumerate() {
                        let wo = &layer.w[o];
                        for (i, p) in prev.iter_mut().enumerate() {
                            *p += wo[i] * d;
                        }
                    }
                    delta = prev;
                }
            }
        }

        // ---- Adam update ----
        self.t += 1;
        let (b1, b2, eps): (f64, f64, f64) = (0.9, 0.999, 1e-8);
        let bc1 = 1.0 - b1.powi(self.t as i32);
        let bc2 = 1.0 - b2.powi(self.t as i32);
        for (li, layer) in self.layers.iter_mut().enumerate() {
            for o in 0..layer.b.len() {
                // weights
                for i in 0..layer.w[o].len() {
                    let g = gw[li][o][i];
                    let m = b1 * layer.mw[o][i] + (1.0 - b1) * g;
                    let v = b2 * layer.vw[o][i] + (1.0 - b2) * g * g;
                    layer.mw[o][i] = m;
                    layer.vw[o][i] = v;
                    let mhat = m / bc1;
                    let vhat = v / bc2;
                    layer.w[o][i] -= lr * mhat / (vhat.sqrt() + eps);
                }
                // bias
                let g = gb[li][o];
                let m = b1 * layer.mb[o] + (1.0 - b1) * g;
                let v = b2 * layer.vb[o] + (1.0 - b2) * g * g;
                layer.mb[o] = m;
                layer.vb[o] = v;
                let mhat = m / bc1;
                let vhat = v / bc2;
                layer.b[o] -= lr * mhat / (vhat.sqrt() + eps);
            }
        }

        mse / n as f64
    }
}

// ---------------------------------------------------------------------------
// The surrogate model: owns the net + the I/O standardisation + a metrics record
// ---------------------------------------------------------------------------

/// A trained surrogate: the `Mlp` plus the **target standardisation**
/// (`mean` / `std`) used to map physical delta to the net's training space and
/// back. Predicts in physical units (metres) by inverting the standardisation.
#[derive(Clone, Debug)]
pub struct SurrogateModel {
    net: Mlp,
    /// Target mean (physical units) used to standardise the training targets.
    y_mean: f64,
    /// Target std-dev (physical units); predictions are `mean + std*net_out`.
    y_std: f64,
    /// Train-set MSE in **standardised** target units (dimensionless).
    pub train_mse: f64,
    /// Held-out test-set MSE in **standardised** target units (dimensionless).
    pub test_mse: f64,
    /// Whether a real training run has populated this model.
    pub trained: bool,
}

impl Default for SurrogateModel {
    fn default() -> Self {
        // An untrained placeholder net so `predict` never panics pre-Train.
        let mut rng = Rng::new(0);
        Self {
            net: Mlp::new(1, &mut rng),
            y_mean: 0.0,
            y_std: 1.0,
            train_mse: f64::NAN,
            test_mse: f64::NAN,
            trained: false,
        }
    }
}

impl SurrogateModel {
    /// Sample the truth on a random design over the input box, split into
    /// train/test, **standardise the target**, and train the MLP. Returns the
    /// fitted model with its train/test MSE recorded.
    ///
    /// - `n_samples` total `(x0,x1)` points (min 8), `test_frac` held out for the
    ///   test MSE;
    /// - `hidden` hidden width, `epochs` full-batch Adam passes, `seed` makes the
    ///   whole run reproducible.
    pub fn train(
        n_samples: usize,
        test_frac: f64,
        hidden: usize,
        epochs: usize,
        seed: u64,
    ) -> Self {
        let n = n_samples.max(8);
        let mut rng = Rng::new(seed);

        // Uniform random design over the [0,1]^2 input box.
        let mut xs: Vec<[f64; 2]> = (0..n).map(|_| [rng.unit(), rng.unit()]).collect();
        // Shuffle (Fisher-Yates) then split.
        for i in (1..xs.len()).rev() {
            let j = (rng.next_u64() % (i as u64 + 1)) as usize;
            xs.swap(i, j);
        }
        let ys_phys: Vec<f64> = xs.iter().map(|x| Truth::deflection(x[0], x[1])).collect();

        let test_frac = test_frac.clamp(0.05, 0.5);
        let n_test = ((n as f64 * test_frac).round() as usize).clamp(1, n - 1);
        let n_train = n - n_test;

        // Standardise the target on the TRAIN split only (no test leakage).
        let mean = ys_phys[..n_train].iter().sum::<f64>() / n_train as f64;
        let var = ys_phys[..n_train]
            .iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>()
            / n_train as f64;
        let std = var.sqrt().max(1e-12);

        let ys_std: Vec<f64> = ys_phys.iter().map(|v| (v - mean) / std).collect();

        let train_x = &xs[..n_train];
        let train_y = &ys_std[..n_train];
        let test_x = &xs[n_train..];
        let test_y = &ys_std[n_train..];

        let mut net = Mlp::new(hidden, &mut rng);
        // A steady learning rate; full-batch Adam converges this tiny problem in
        // a few hundred epochs. Epoch count is user/agent-controllable.
        let lr = 0.01;
        let epochs = epochs.clamp(1, 20_000);
        for _ in 0..epochs {
            net.train_step(train_x, train_y, lr);
        }

        // Final metrics (the post-training values).
        let mse = |net: &Mlp, xs: &[[f64; 2]], ys: &[f64]| -> f64 {
            if xs.is_empty() {
                return f64::NAN;
            }
            xs.iter()
                .zip(ys)
                .map(|(x, &y)| {
                    let e = net.predict(x) - y;
                    e * e
                })
                .sum::<f64>()
                / xs.len() as f64
        };
        let train_mse_final = mse(&net, train_x, train_y);
        let test_mse_final = mse(&net, test_x, test_y);

        Self {
            net,
            y_mean: mean,
            y_std: std,
            train_mse: train_mse_final,
            test_mse: test_mse_final,
            trained: true,
        }
    }

    /// Predict the surrogate's deflection `delta` (physical metres) at normalised
    /// inputs `(x0, x1)`. Microsecond-cheap — this is what drives the real-time
    /// slider readout. Inverts the target standardisation.
    pub fn predict(&self, x0: f64, x1: f64) -> f64 {
        let out = self.net.predict(&[x0.clamp(0.0, 1.0), x1.clamp(0.0, 1.0)]);
        self.y_mean + self.y_std * out
    }
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the Surrogate Model workbench: the training hyper-
/// parameters, the two live prediction inputs, and the most recently trained
/// surrogate.
pub struct SurrogateWorkbenchState {
    /// Number of `(x0,x1)` points sampled from the truth for training+test.
    pub n_samples: usize,
    /// Fraction of samples held out for the test-MSE (`0.05 .. 0.5`).
    pub test_frac: f64,
    /// Hidden-layer width of the `2->H->H->1` MLP.
    pub hidden: usize,
    /// Full-batch Adam training passes.
    pub epochs: usize,
    /// Live prediction input #0 — the normalised **load** factor `x0` in `[0,1]`.
    pub input0: f64,
    /// Live prediction input #1 — the normalised **length** factor `x1` in
    /// `[0,1]`.
    pub input1: f64,
    /// The most recently trained surrogate (untrained until **Train** runs).
    pub model: SurrogateModel,
}

impl Default for SurrogateWorkbenchState {
    fn default() -> Self {
        let mut s = Self {
            n_samples: 400,
            test_frac: 0.2,
            hidden: 16,
            epochs: 600,
            input0: 0.5,
            input1: 0.5,
            model: SurrogateModel::default(),
        };
        // Train once on construction so the panel shows a working surrogate
        // immediately (and the readout has metrics from the first frame).
        s.train();
        s
    }
}

impl SurrogateWorkbenchState {
    /// Sample the truth and (re)train the surrogate with the current hyper-
    /// parameters. Shared by the in-panel **Train** button and the
    /// `surrogate.train` bridge id so both run the SAME path. A fixed seed makes
    /// each Train reproducible.
    pub fn train(&mut self) {
        self.model = SurrogateModel::train(
            self.n_samples,
            self.test_frac,
            self.hidden,
            self.epochs,
            0xC0FF_EE15_F00D,
        );
    }

    /// The true (closed-form) deflection at the current input sliders.
    pub fn true_value(&self) -> f64 {
        Truth::deflection(self.input0, self.input1)
    }

    /// The surrogate's predicted deflection at the current input sliders.
    pub fn predicted_value(&self) -> f64 {
        self.model.predict(self.input0, self.input1)
    }

    /// Captions of every control the agent bridge can `SetControl`.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "Load factor",
            "Length factor",
            "Training samples",
            "Hidden width",
            "Epochs",
            "Test fraction",
        ]
    }

    /// Set one labelled control by caption for the agent `SetControl` bridge.
    /// Fail-loud on an unknown caption / wrong type / out-of-range value; no
    /// state is written on error and nothing panics. All controls are finite
    /// numbers; the count controls are rounded to a sensible integer.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        let num = |value: &crate::agent_commands::AgentValue| -> Result<f64, String> {
            let v = value.as_f64()?;
            if !v.is_finite() {
                return Err(format!("{name}: value must be finite, got {v}"));
            }
            Ok(v)
        };
        match name {
            "Load factor" => {
                self.input0 = num(value)?.clamp(0.0, 1.0);
                Ok(())
            }
            "Length factor" => {
                self.input1 = num(value)?.clamp(0.0, 1.0);
                Ok(())
            }
            "Training samples" => {
                let v = num(value)?;
                if v < 8.0 {
                    return Err(format!("Training samples: must be >= 8, got {v}"));
                }
                self.n_samples = v.round() as usize;
                Ok(())
            }
            "Hidden width" => {
                let v = num(value)?;
                if v < 1.0 {
                    return Err(format!("Hidden width: must be >= 1, got {v}"));
                }
                self.hidden = (v.round() as usize).min(256);
                Ok(())
            }
            "Epochs" => {
                let v = num(value)?;
                if v < 1.0 {
                    return Err(format!("Epochs: must be >= 1, got {v}"));
                }
                self.epochs = (v.round() as usize).min(20_000);
                Ok(())
            }
            "Test fraction" => {
                let v = num(value)?;
                if !(0.05..=0.5).contains(&v) {
                    return Err(format!("Test fraction: must be in [0.05, 0.5], got {v}"));
                }
                self.test_frac = v;
                Ok(())
            }
            other => Err(format!("unknown surrogate control: {other:?}")),
        }
    }

    /// Readout for the agent `ReadReadout` bridge: train/test MSE and the live
    /// surrogate-vs-true prediction at the input sliders. Always `Some`.
    pub fn agent_readout(&self) -> Option<String> {
        let m = &self.model;
        let metrics = if m.trained {
            format!(
                "train MSE={:.3e}, test MSE={:.3e} (standardised target units)",
                m.train_mse, m.test_mse
            )
        } else {
            "(not trained)".to_string()
        };
        let truth = self.true_value();
        let pred = self.predicted_value();
        let rel = if truth.abs() > 1e-15 {
            (pred - truth).abs() / truth.abs() * 100.0
        } else {
            0.0
        };
        Some(format!(
            "Surrogate (cantilever \u{03B4}=P\u{00B7}L\u{00B3}/3EI) \u{00B7} {metrics} \u{00B7} \
             @(load={:.3}, length={:.3}): true={:.6e} m, surrogate={:.6e} m (err {:.2}%)",
            self.input0, self.input1, truth, pred, rel
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (train)
// ---------------------------------------------------------------------------

/// Run the sample-then-train (the in-panel **Train** action). Factored out so
/// the button and the `surrogate.train` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.surrogate.train();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Surrogate Model workbench. A no-op unless toggled on via View ->
/// Surrogate Model.
pub fn draw_surrogate_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_surrogate_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_surrogate_workbench",
        "Surrogate Model (ML emulator for a solver)",
        surrogate_workbench_body,
    );
    if close {
        app.show_surrogate_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn surrogate_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house ML surrogate [a surrogate model trains a small neural net on samples from an \
             expensive solver, then predicts the output INSTANTLY when inputs change \u{2014} the \
             slider-moves-and-result-updates-instantly idea. Here the ground truth is the closed-form \
             cantilever tip deflection \u{03B4} = P\u{00B7}L\u{00B3}/(3\u{00B7}E\u{00B7}I) over two \
             normalised inputs (load, length), so the surrogate's accuracy is validatable against \
             truth. Set the hyper-parameters, press Train (samples the truth + fits a \
             2\u{2192}H\u{2192}H\u{2192}1 MLP with ReLU + Adam on MSE), then move the input sliders: \
             the surrogate prediction is shown live next to the true value].",
        )
        .weak()
        .small(),
    );
    ui.separator();

    // --- Training hyper-parameters + Train button ---------------------------
    ui.label(egui::RichText::new("Training").strong());
    egui::Grid::new("surrogate_train_params")
        .num_columns(2)
        .spacing([12.0, 4.0])
        .show(ui, |ui| {
            int_row(
                ui,
                "Training samples",
                "number of (load, length) points sampled from the truth",
                &mut app.surrogate.n_samples,
                8,
                4000,
                5.0,
            );
            int_row(
                ui,
                "Hidden width",
                "width H of the 2\u{2192}H\u{2192}H\u{2192}1 MLP hidden layers",
                &mut app.surrogate.hidden,
                1,
                256,
                1.0,
            );
            int_row(
                ui,
                "Epochs",
                "full-batch Adam training passes",
                &mut app.surrogate.epochs,
                1,
                20_000,
                10.0,
            );
            // Test fraction (a real f64 control).
            let lbl = ui.label("Test fraction");
            ui.add(
                egui::DragValue::new(&mut app.surrogate.test_frac)
                    .speed(0.01)
                    .range(0.05..=0.5)
                    .max_decimals(3),
            )
            .labelled_by(lbl.id)
            .on_hover_text("fraction of samples held out for the test MSE");
            ui.end_row();
        });

    ui.add_space(4.0);
    if ui
        .button("\u{25B6} Train surrogate")
        .on_hover_text(
            "Sample the true cantilever solver and fit the MLP (ReLU + Adam on MSE). Reports \
             train / test MSE.",
        )
        .clicked()
    {
        app.surrogate.train();
    }

    // --- Training metrics ---------------------------------------------------
    let m = &app.surrogate.model;
    if m.trained {
        ui.label(
            egui::RichText::new(format!(
                "Train MSE = {:.3e}   \u{00B7}   Test MSE = {:.3e}   (standardised target units)",
                m.train_mse, m.test_mse
            ))
            .strong(),
        );
        ui.label(
            egui::RichText::new(
                "(MSE is in zero-mean/unit-variance target space, so < 1 means the net beats a \
                 constant-mean predictor; a small value means it learned the function.)",
            )
            .weak()
            .small(),
        );
    } else {
        ui.label(egui::RichText::new("Press Train to fit the surrogate.").weak());
    }

    ui.separator();

    // --- Real-time prediction: sliders -> instant surrogate vs true ---------
    ui.label(egui::RichText::new("Real-time prediction (what-if)").strong());
    egui::Grid::new("surrogate_predict_inputs")
        .num_columns(2)
        .spacing([12.0, 4.0])
        .show(ui, |ui| {
            let lbl0 = ui.label("Load factor");
            ui.add(
                egui::Slider::new(&mut app.surrogate.input0, 0.0..=1.0)
                    .text("x\u{2080}")
                    .max_decimals(3),
            )
            .labelled_by(lbl0.id)
            .on_hover_text(format!(
                "normalised load \u{2192} P = {:.0} N",
                Truth::load(app.surrogate.input0)
            ));
            ui.end_row();

            let lbl1 = ui.label("Length factor");
            ui.add(
                egui::Slider::new(&mut app.surrogate.input1, 0.0..=1.0)
                    .text("x\u{2081}")
                    .max_decimals(3),
            )
            .labelled_by(lbl1.id)
            .on_hover_text(format!(
                "normalised length \u{2192} L = {:.3} m",
                Truth::length(app.surrogate.input1)
            ));
            ui.end_row();
        });

    // The live surrogate-vs-true readout (this is the instant what-if result).
    let truth = app.surrogate.true_value();
    let pred = app.surrogate.predicted_value();
    let rel = if truth.abs() > 1e-15 {
        (pred - truth).abs() / truth.abs() * 100.0
    } else {
        0.0
    };
    ui.add_space(4.0);
    egui::Grid::new("surrogate_predict_out")
        .num_columns(2)
        .spacing([12.0, 2.0])
        .show(ui, |ui| {
            ui.label("True deflection \u{03B4}");
            ui.label(
                egui::RichText::new(format!("{truth:.6e} m"))
                    .monospace()
                    .color(egui::Color32::from_rgb(150, 200, 150)),
            );
            ui.end_row();
            ui.label("Surrogate prediction");
            ui.label(
                egui::RichText::new(format!("{pred:.6e} m"))
                    .monospace()
                    .color(egui::Color32::from_rgb(150, 190, 240)),
            );
            ui.end_row();
            ui.label("Relative error");
            ui.label(egui::RichText::new(format!("{rel:.2} %")).monospace());
            ui.end_row();
        });

    // --- Sweep plot: surrogate vs true over the length input ----------------
    ui.add_space(6.0);
    ui.label(
        egui::RichText::new("Surrogate vs true sweep (vary length at the current load)").strong(),
    );
    let x0 = app.surrogate.input0;
    let model = &app.surrogate.model;
    let true_pts: PlotPoints = (0..=120)
        .map(|i| {
            let x1 = i as f64 / 120.0;
            [Truth::length(x1), Truth::deflection(x0, x1)]
        })
        .collect();
    let surr_pts: PlotPoints = (0..=120)
        .map(|i| {
            let x1 = i as f64 / 120.0;
            [Truth::length(x1), model.predict(x0, x1)]
        })
        .collect();
    Plot::new("surrogate_sweep")
        .height(200.0)
        .legend(Legend::default())
        .show(ui, |pui| {
            pui.line(Line::new(true_pts).name("true \u{03B4}(L)"));
            pui.line(Line::new(surr_pts).name("surrogate"));
        });
}

/// One labelled integer `DragValue` parameter row inside a grid. The caption is
/// a named label the DragValue is `labelled_by`, so the agent bridge / a screen
/// reader can find the spin button by its caption text (the AI-drivable name).
fn int_row(
    ui: &mut egui::Ui,
    caption: &str,
    hover: &str,
    value: &mut usize,
    lo: usize,
    hi: usize,
    speed: f64,
) {
    let lbl = ui.label(caption);
    ui.add(
        egui::DragValue::new(value)
            .speed(speed)
            .range(lo..=hi)
            .max_decimals(0),
    )
    .labelled_by(lbl.id)
    .on_hover_text(hover);
    ui.end_row();
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    #[test]
    fn truth_matches_closed_form() {
        // Spot-check delta = P L^3 / (3 E I) against a hand computation.
        let x0 = 1.0; // P = 1100 N
        let x1 = 1.0; // L = 2.5 m
        let p = 1100.0;
        let l = 2.5;
        let expected = p * l * l * l / (3.0 * 200e9 * 8.0e-6);
        assert!((Truth::deflection(x0, x1) - expected).abs() < 1e-15);
        // Zero-input maps to the range minimums (P=100, L=0.5), nonzero delta.
        assert!(Truth::deflection(0.0, 0.0) > 0.0);
        // Monotonic in both inputs (more load / longer beam -> more deflection).
        assert!(Truth::deflection(0.8, 0.5) > Truth::deflection(0.2, 0.5));
        assert!(Truth::deflection(0.5, 0.8) > Truth::deflection(0.5, 0.2));
    }

    #[test]
    fn surrogate_learns_the_function_small_test_mse() {
        // THE validation: after training, the held-out test MSE (standardised
        // target units) must be small — the net actually learned delta(x0,x1).
        let m = SurrogateModel::train(500, 0.2, 16, 800, 1234);
        assert!(m.trained);
        assert!(
            m.train_mse.is_finite() && m.test_mse.is_finite(),
            "metrics must be finite: train={}, test={}",
            m.train_mse,
            m.test_mse
        );
        // In zero-mean/unit-variance target space a constant-mean predictor
        // scores ~1.0; a net that learned the surface must beat it by a wide
        // margin. 0.02 is a comfortable, non-flaky bound for this smooth 2-input
        // function with this budget.
        assert!(
            m.test_mse < 0.02,
            "surrogate test MSE should be small (learned the function); got {}",
            m.test_mse
        );
        // Train MSE should also be small.
        assert!(
            m.train_mse < 0.02,
            "surrogate train MSE should be small; got {}",
            m.train_mse
        );
    }

    #[test]
    fn surrogate_prediction_tracks_truth_in_physical_units() {
        // After training, the physical-unit prediction must be close to truth at
        // arbitrary query points (relative error well under 10%).
        let m = SurrogateModel::train(600, 0.2, 24, 1000, 7);
        for &(x0, x1) in &[(0.2, 0.3), (0.5, 0.5), (0.85, 0.7), (0.1, 0.95)] {
            let truth = Truth::deflection(x0, x1);
            let pred = m.predict(x0, x1);
            let rel = (pred - truth).abs() / truth.abs();
            assert!(
                rel < 0.10,
                "surrogate within 10% at ({x0},{x1}): true={truth:.4e}, pred={pred:.4e}, rel={rel:.3}"
            );
        }
    }

    #[test]
    fn training_is_deterministic_for_a_seed() {
        let a = SurrogateModel::train(200, 0.2, 8, 100, 42);
        let b = SurrogateModel::train(200, 0.2, 8, 100, 42);
        assert_eq!(a.test_mse.to_bits(), b.test_mse.to_bits());
        assert_eq!(
            a.predict(0.33, 0.66).to_bits(),
            b.predict(0.33, 0.66).to_bits()
        );
    }

    #[test]
    fn predict_is_fast_enough_for_realtime() {
        // Sanity: a single predict is microsecond-scale (the whole point — it can
        // run every frame as a slider moves). 10k predicts must be near-instant.
        let m = SurrogateModel::train(200, 0.2, 16, 200, 1);
        let t0 = std::time::Instant::now();
        let mut acc = 0.0;
        for i in 0..10_000 {
            let x = i as f64 / 10_000.0;
            acc += m.predict(x, 1.0 - x);
        }
        let dt = t0.elapsed();
        assert!(acc.is_finite());
        assert!(
            dt.as_millis() < 500,
            "10k surrogate predicts should be well under 500ms; took {dt:?}"
        );
    }

    #[test]
    fn default_state_is_trained() {
        let s = SurrogateWorkbenchState::default();
        assert!(s.model.trained, "default state trains a surrogate eagerly");
        assert!(s.model.test_mse.is_finite());
    }

    #[test]
    fn agent_set_inputs_and_hyperparams_round_trip() {
        let mut s = SurrogateWorkbenchState::default();
        s.agent_set("Load factor", &AgentValue::Float(0.25))
            .unwrap();
        assert_eq!(s.input0, 0.25);
        s.agent_set("Length factor", &AgentValue::Float(0.75))
            .unwrap();
        assert_eq!(s.input1, 0.75);
        // Out-of-range inputs clamp to [0,1] (no error — they are sliders).
        s.agent_set("Load factor", &AgentValue::Float(2.0)).unwrap();
        assert_eq!(s.input0, 1.0);
        s.agent_set("Hidden width", &AgentValue::Float(32.0))
            .unwrap();
        assert_eq!(s.hidden, 32);
        s.agent_set("Epochs", &AgentValue::Float(300.0)).unwrap();
        assert_eq!(s.epochs, 300);
        s.agent_set("Training samples", &AgentValue::Float(256.0))
            .unwrap();
        assert_eq!(s.n_samples, 256);
        s.agent_set("Test fraction", &AgentValue::Float(0.3))
            .unwrap();
        assert_eq!(s.test_frac, 0.3);
        // Fail-loud cases.
        assert!(s.agent_set("Epochs", &AgentValue::Float(0.0)).is_err());
        assert!(s
            .agent_set("Test fraction", &AgentValue::Float(0.9))
            .is_err());
        assert!(s
            .agent_set("Training samples", &AgentValue::Float(4.0))
            .is_err());
        assert!(s.agent_set("Mass", &AgentValue::Float(1.0)).is_err());
        assert!(s
            .agent_set("Hidden width", &AgentValue::Float(f64::NAN))
            .is_err());
    }

    #[test]
    fn control_names_listed_and_nonempty() {
        let names = SurrogateWorkbenchState::agent_control_names();
        assert!(names.contains(&"Load factor"));
        assert!(names.contains(&"Length factor"));
        assert!(names.contains(&"Epochs"));
        assert!(names.contains(&"Hidden width"));
    }

    #[test]
    fn readout_reports_mse_and_prediction() {
        let s = SurrogateWorkbenchState::default();
        let r = s.agent_readout().expect("readout always present");
        assert!(r.contains("train MSE="), "got: {r}");
        assert!(r.contains("test MSE="), "got: {r}");
        assert!(r.contains("surrogate="), "got: {r}");
        assert!(r.contains("true="), "got: {r}");
    }

    #[test]
    fn run_bridge_helper_trains_through_app() {
        let mut app = ValenxApp::default();
        app.surrogate
            .agent_set("Epochs", &AgentValue::Float(120.0))
            .unwrap();
        app.surrogate
            .agent_set("Training samples", &AgentValue::Float(150.0))
            .unwrap();
        run(&mut app);
        assert!(app.surrogate.model.trained);
        let r = app.surrogate.agent_readout().unwrap();
        assert!(r.contains("test MSE="), "got: {r}");
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
            draw_surrogate_workbench(app, ctx);
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
        assert!(!app.show_surrogate_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_surrogate_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_surrogate_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every parameter DragValue/Slider is a SpinButton/Slider and must be
        // `labelled_by` its caption so an AI / screen reader can find it by
        // caption text (the AI-drivable name).
        let mut app = ValenxApp::default();
        app.show_surrogate_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let controls: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| matches!(n.role(), Role::SpinButton | Role::Slider))
            .collect();
        assert!(
            !controls.is_empty(),
            "expected the parameter numeric controls as spin buttons / sliders"
        );
        assert!(
            controls.iter().all(|n| {
                !n.labelled_by().is_empty() || n.name().is_some_and(|s| !s.trim().is_empty())
            }),
            "every numeric control must be named or labelled_by a caption (AI-drivable)"
        );
        assert!(
            controls.iter().all(|n| {
                n.name().is_some_and(|s| !s.trim().is_empty())
                    || n.labelled_by()
                        .iter()
                        .any(|id| by_id.get(id).is_some_and(|t| t.name().is_some()))
            }),
            "every numeric control resolves to a named caption node"
        );
        // The live-input captions are named nodes.
        assert!(
            has_named_node(&nodes, "Load factor"),
            "'Load factor' caption is a named node"
        );
        assert!(
            has_named_node(&nodes, "Length factor"),
            "'Length factor' caption is a named node"
        );
        assert!(
            has_named_node(&nodes, "Epochs"),
            "'Epochs' caption is a named node"
        );
    }
}
