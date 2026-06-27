//! The right-side **Bond Graph** workbench — a native, in-house **multi-domain
//! systems-modelling** tool. A *bond graph* describes a physical system by how
//! **power** flows through it: every bond carries an **effort** `e` and a
//! **flow** `f` whose product `e·f` is power, and the SAME small element set
//! models mechanical, electrical, hydraulic and thermal sub-systems uniformly
//! (force·velocity, voltage·current, pressure·flow, …). That domain-independence
//! is the whole point: a mass, an inductor and a fluid inertance are all the one
//! element **I** (inertance) on their respective bonds.
//!
//! ## What V1 ships (preset-based, with a real graph→state derivation)
//!
//! The standard bond-graph **elements** are exposed as canvas node kinds and
//! drawn on the SAME in-house node-graph canvas as [`crate::nodegraph_workbench`]
//! (drag nodes, see the wired power bonds), so the bond graph is *visualised*:
//! the 1-ports **R** (resistance / damper), **C** (compliance / spring /
//! capacitor), **I** (inertance / mass / inductor), **Se** (effort source) and
//! **Sf** (flow source); the 2-ports **TF** (transformer) and **GY** (gyrator);
//! and the **0-junction** (common effort) and **1-junction** (common flow).
//!
//! For the **state equations** V1 takes the well-trodden *preset* route: three
//! canonical systems (mass–spring–damper, series RLC, DC motor) ship with their
//! bond graph AND the **linear state-space `dx/dt = A·x + B·u`** that the
//! standard bond-graph derivation yields, where the state variables are exactly
//! the **energy variables stored on the I and C elements** — generalised
//! momentum `p` on each **I**, generalised displacement `q` on each **C** — under
//! **integral causality** (the bond-graph normal form). Those `A`/`B` are derived
//! by hand from the junction equations (see [`BondGraphPreset::state_space`]) and
//! are unit-tested against the systems' independently-known analytic ODEs
//! (mass–spring–damper → `m x'' + b x' + k x = F`; series RLC →
//! `L q'' + R q' + q/C = V`). The arbitrary-graph causality-assignment +
//! symbolic state-equation derivation is a deeper follow-up; this is the honest
//! *"preset-based V1, general derivation later"* scope.
//!
//! The state ODEs are integrated with classical **RK4** and the response is
//! plotted (egui_plot).
//!
//! ## AI-drivable surface (the #1 standing release gate)
//!
//! Mirrors every other workbench: a [`crate::workbench_chrome::workbench_shell`]
//! panel gated on [`crate::ValenxApp::show_bondgraph_workbench`], toggled from
//! the View menu and openable by the agent bridge under the id `"bondgraph"`.
//! The bridge can:
//! - pick the preset and set every element parameter through labelled controls
//!   (`agent_set` / `agent_control_names`): mass / damping / stiffness / force,
//!   resistance / inductance / capacitance / voltage, the motor constants, plus
//!   the simulation duration;
//! - read a status line (`agent_readout`): the derived ODE order, the natural
//!   frequency and damping ratio, and the final state;
//! - fire the solve via the `RunCommand` id `bondgraph.solve`, which runs the
//!   SAME derive-then-integrate path the in-panel **Solve** button calls.

use eframe::egui;
use egui_plot::{Legend, Line, PlotPoints};

use crate::plot_ui::managed_plot_mem_cfg;
use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Bond-graph elements (the canvas node kinds — visualisation layer)
// ---------------------------------------------------------------------------

/// A standard bond-graph element. These are the node kinds drawn on the canvas
/// so the user *sees* the bond graph of the selected preset; the numerical
/// state-space lives in [`BondGraphPreset`]. Power-bond semantics (effort `e`,
/// flow `f`, power `e·f`) are domain-independent — the same element models a
/// mechanical, electrical, hydraulic or thermal port.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BgElement {
    /// **R** — resistance / damper. Dissipates power: `e = R·f` (electrical
    /// resistor, mechanical damper, fluid restriction).
    R,
    /// **C** — compliance. Stores displacement energy: `e = q/C`, `q = ∫f dt`
    /// (spring `C = 1/k`, electrical capacitor, fluid compliance).
    C,
    /// **I** — inertance. Stores momentum energy: `f = p/I`, `p = ∫e dt`
    /// (mass / inductor / fluid inertance).
    I,
    /// **Se** — effort source (force / voltage / pressure source).
    Se,
    /// **Sf** — flow source (velocity / current / volumetric-flow source).
    Sf,
    /// **TF** — transformer. Power-conserving effort/flow scaling by a modulus
    /// `m`: `e₁ = m·e₂`, `f₂ = m·f₁` (lever, gear, ideal electrical transformer).
    TF,
    /// **GY** — gyrator. Power-conserving effort↔flow cross-coupling by `r`:
    /// `e₁ = r·f₂`, `e₂ = r·f₁` (the DC-motor back-EMF / torque coupling).
    GY,
    /// **0-junction** — common **effort** (all bonds share `e`; flows sum to
    /// zero). The electrical *node* / mechanical *equal-force* connection.
    J0,
    /// **1-junction** — common **flow** (all bonds share `f`; efforts sum to
    /// zero). The mechanical *equal-velocity* / electrical *series-loop*
    /// connection.
    J1,
}

impl BgElement {
    /// Every element, in palette order.
    pub const ALL: [BgElement; 9] = [
        BgElement::R,
        BgElement::C,
        BgElement::I,
        BgElement::Se,
        BgElement::Sf,
        BgElement::TF,
        BgElement::GY,
        BgElement::J0,
        BgElement::J1,
    ];

    /// Short symbol drawn on the node body (`R`, `C`, `I`, `Se`, `0`, `1`, …).
    pub fn symbol(self) -> &'static str {
        match self {
            BgElement::R => "R",
            BgElement::C => "C",
            BgElement::I => "I",
            BgElement::Se => "Se",
            BgElement::Sf => "Sf",
            BgElement::TF => "TF",
            BgElement::GY => "GY",
            BgElement::J0 => "0",
            BgElement::J1 => "1",
        }
    }

    /// Longer human label for hovers / legends.
    pub fn label(self) -> &'static str {
        match self {
            BgElement::R => "R (resistance / damper)",
            BgElement::C => "C (compliance / spring / capacitor)",
            BgElement::I => "I (inertance / mass / inductor)",
            BgElement::Se => "Se (effort source)",
            BgElement::Sf => "Sf (flow source)",
            BgElement::TF => "TF (transformer)",
            BgElement::GY => "GY (gyrator)",
            BgElement::J0 => "0-junction (common effort)",
            BgElement::J1 => "1-junction (common flow)",
        }
    }

    /// Colour tint for the node header (visual grouping by element class).
    fn color(self) -> egui::Color32 {
        match self {
            BgElement::R => egui::Color32::from_rgb(120, 80, 70),
            BgElement::C => egui::Color32::from_rgb(70, 110, 90),
            BgElement::I => egui::Color32::from_rgb(70, 95, 130),
            BgElement::Se | BgElement::Sf => egui::Color32::from_rgb(120, 100, 60),
            BgElement::TF | BgElement::GY => egui::Color32::from_rgb(100, 75, 120),
            BgElement::J0 | BgElement::J1 => egui::Color32::from_rgb(70, 78, 95),
        }
    }
}

/// One placed bond-graph node: an element, a label (e.g. `"I: m"`) and a canvas
/// position. Pure data so the preset graphs are trivially constructible /
/// testable headless.
#[derive(Clone, Debug)]
pub struct BgNode {
    /// The element this node represents.
    pub element: BgElement,
    /// Short caption shown under the symbol (e.g. `"m"`, `"b"`, `"1/k"`).
    pub caption: String,
    /// Canvas-local top-left position.
    pub pos: egui::Pos2,
}

/// A power **bond** between two placed nodes (drawn as a wire on the canvas).
/// Bonds are undirected for drawing purposes here (the half-arrow / causal
/// stroke is a later refinement); index `from`/`to` into the preset node list.
#[derive(Clone, Copy, Debug)]
pub struct BgBond {
    /// Index of the source node in the preset node list.
    pub from: usize,
    /// Index of the destination node.
    pub to: usize,
}

// ---------------------------------------------------------------------------
// Presets — the bond graph + its derived linear state-space
// ---------------------------------------------------------------------------

/// The three canonical preset systems V1 ships. Each owns its element
/// parameters, builds its **bond graph** (for the canvas) and derives its
/// **linear state-space** `dx/dt = A·x + B·u` from the junction equations under
/// integral causality.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondGraphPreset {
    /// Mechanical mass–spring–damper: `Se:F — 1 — {I:m, R:b, C:1/k}`.
    MassSpringDamper,
    /// Electrical series RLC: `Se:V — 1 — {I:L, R:R, C:C}`.
    Rlc,
    /// DC motor: `Se:V — 1ₑ — {R:Rₐ, I:Lₐ} — GY:Km — 1ₘ — {I:J, R:b}`.
    DcMotor,
}

impl BondGraphPreset {
    /// Every preset, in selector order.
    pub const ALL: [BondGraphPreset; 3] = [
        BondGraphPreset::MassSpringDamper,
        BondGraphPreset::Rlc,
        BondGraphPreset::DcMotor,
    ];

    /// Selector / menu label.
    pub fn label(self) -> &'static str {
        match self {
            BondGraphPreset::MassSpringDamper => "Mass-spring-damper",
            BondGraphPreset::Rlc => "Series RLC circuit",
            BondGraphPreset::DcMotor => "DC motor",
        }
    }

    /// Stable id used by the agent bridge.
    pub fn id(self) -> &'static str {
        match self {
            BondGraphPreset::MassSpringDamper => "msd",
            BondGraphPreset::Rlc => "rlc",
            BondGraphPreset::DcMotor => "dcmotor",
        }
    }

    /// Parse a preset from an id / friendly alias (case-insensitive).
    pub fn from_id(s: &str) -> Option<BondGraphPreset> {
        match s.trim().to_ascii_lowercase().as_str() {
            "msd" | "massspringdamper" | "mass-spring-damper" | "spring" | "mechanical" => {
                Some(BondGraphPreset::MassSpringDamper)
            }
            "rlc" | "circuit" | "electrical" => Some(BondGraphPreset::Rlc),
            "dcmotor" | "motor" | "dc" => Some(BondGraphPreset::DcMotor),
            _ => None,
        }
    }
}

/// All editable element parameters for every preset (only the active preset's
/// fields are read by `Self::state_space` / the UI, but holding them all lets
/// the user flip presets without losing values). SI units throughout.
#[derive(Clone, Copy, Debug)]
pub struct BondGraphParams {
    // --- mass-spring-damper ---
    /// Mass `m` (kg) — the I element.
    pub mass: f64,
    /// Damping `b` (N·s/m) — the R element.
    pub damping: f64,
    /// Stiffness `k` (N/m) — the C element is the compliance `1/k`.
    pub stiffness: f64,
    /// Applied force `F` (N) — the Se source (a step input).
    pub force: f64,

    // --- series RLC ---
    /// Resistance `R` (Ω).
    pub resistance: f64,
    /// Inductance `L` (H) — the I element.
    pub inductance: f64,
    /// Capacitance `C` (F) — the C element.
    pub capacitance: f64,
    /// Source voltage `V` (V) — the Se source (a step input).
    pub voltage: f64,

    // --- DC motor ---
    /// Armature resistance `Rₐ` (Ω).
    pub r_arm: f64,
    /// Armature inductance `Lₐ` (H).
    pub l_arm: f64,
    /// Rotor inertia `J` (kg·m²).
    pub inertia: f64,
    /// Rotor viscous friction `b` (N·m·s).
    pub friction: f64,
    /// Motor / back-EMF constant `Km` (N·m/A = V·s/rad) — the GY modulus.
    pub k_motor: f64,
    /// Supply voltage `Vₘ` (V) — the Se source (a step input).
    pub motor_voltage: f64,

    // --- simulation ---
    /// Total simulated time (s).
    pub duration: f64,
}

impl Default for BondGraphParams {
    fn default() -> Self {
        Self {
            // Mass-spring-damper: m=1, b=0.5, k=20 → ωn≈4.47 rad/s, ζ≈0.056
            // (lightly damped — a clear ringing step response).
            mass: 1.0,
            damping: 0.5,
            stiffness: 20.0,
            force: 1.0,
            // Series RLC: R=2, L=1, C=0.01 → ωn=10 rad/s, ζ=0.1.
            resistance: 2.0,
            inductance: 1.0,
            capacitance: 0.01,
            voltage: 1.0,
            // DC motor: typical small-motor scales.
            r_arm: 1.0,
            l_arm: 0.5,
            inertia: 0.01,
            friction: 0.1,
            k_motor: 0.1,
            motor_voltage: 12.0,
            duration: 5.0,
        }
    }
}

/// A linear time-invariant state-space realisation `dx/dt = A·x + B·u` with a
/// constant (step) scalar input `u`. The dimension is the **ODE order** = the
/// number of independent energy stores (I + C elements) under integral
/// causality. Kept tiny (≤2 states for V1's presets) and `Copy`-free `Vec`s so
/// the integrator and the analytic-ODE tests share one type.
#[derive(Clone, Debug)]
pub struct StateSpace {
    /// `n×n` system matrix, row-major.
    pub a: Vec<Vec<f64>>,
    /// `n` input column (multiplies the scalar step input `u`).
    pub b: Vec<f64>,
    /// The constant step input value `u`.
    pub u: f64,
    /// Human captions for each state variable (for the plot legend / readout).
    pub state_names: Vec<&'static str>,
}

impl StateSpace {
    /// ODE order = number of state variables = number of I/C energy stores.
    pub fn order(&self) -> usize {
        self.b.len()
    }

    /// Evaluate `dx/dt = A·x + B·u` at state `x`.
    fn deriv(&self, x: &[f64]) -> Vec<f64> {
        let n = self.b.len();
        let mut dx = vec![0.0; n];
        for (i, dxi) in dx.iter_mut().enumerate() {
            let mut acc = self.b[i] * self.u;
            for (j, &xj) in x.iter().enumerate().take(n) {
                acc += self.a[i][j] * xj;
            }
            *dxi = acc;
        }
        dx
    }

    /// Integrate from rest (`x = 0`) over `[0, t_end]` with classical **RK4** at
    /// fixed step `dt`, returning `(t, x)` samples (including `t = 0`). A
    /// non-finite or non-positive `dt`/`t_end` yields just the initial sample
    /// (fail-soft, never a hang).
    pub fn integrate_rk4(&self, t_end: f64, dt: f64) -> Vec<(f64, Vec<f64>)> {
        let n = self.b.len();
        let mut out = vec![(0.0, vec![0.0; n])];
        if !(dt.is_finite() && t_end.is_finite()) || dt <= 0.0 || t_end <= 0.0 {
            return out;
        }
        let steps = (t_end / dt).ceil() as usize;
        let steps = steps.min(2_000_000); // hard cap — never spin forever
        let mut x = vec![0.0; n];
        let mut t = 0.0;
        let add = |a: &[f64], b: &[f64], s: f64| -> Vec<f64> {
            a.iter().zip(b).map(|(ai, bi)| ai + s * bi).collect()
        };
        for _ in 0..steps {
            let k1 = self.deriv(&x);
            let k2 = self.deriv(&add(&x, &k1, dt * 0.5));
            let k3 = self.deriv(&add(&x, &k2, dt * 0.5));
            let k4 = self.deriv(&add(&x, &k3, dt));
            for i in 0..n {
                x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            t += dt;
            out.push((t, x.clone()));
        }
        out
    }
}

impl BondGraphPreset {
    /// Build the **bond graph** (nodes + power bonds) for this preset's canvas
    /// drawing, captioned with the current parameter symbols. Pure data.
    pub fn bond_graph(self, p: &BondGraphParams) -> (Vec<BgNode>, Vec<BgBond>) {
        let node = |element: BgElement, caption: &str, x: f32, y: f32| BgNode {
            element,
            caption: caption.to_string(),
            pos: egui::pos2(x, y),
        };
        match self {
            BondGraphPreset::MassSpringDamper => {
                // Se:F — 1 — {I:m, R:b, C:1/k}
                let nodes = vec![
                    node(BgElement::Se, &format!("F={:.3}", p.force), 20.0, 110.0),
                    node(BgElement::J1, "v", 160.0, 110.0),
                    node(BgElement::I, &format!("m={:.3}", p.mass), 300.0, 30.0),
                    node(BgElement::R, &format!("b={:.3}", p.damping), 300.0, 110.0),
                    node(
                        BgElement::C,
                        &format!("1/k={:.4}", 1.0 / p.stiffness.max(1e-12)),
                        300.0,
                        190.0,
                    ),
                ];
                let bonds = vec![
                    BgBond { from: 0, to: 1 },
                    BgBond { from: 1, to: 2 },
                    BgBond { from: 1, to: 3 },
                    BgBond { from: 1, to: 4 },
                ];
                (nodes, bonds)
            }
            BondGraphPreset::Rlc => {
                // Se:V — 1 — {I:L, R:R, C:C}
                let nodes = vec![
                    node(BgElement::Se, &format!("V={:.3}", p.voltage), 20.0, 110.0),
                    node(BgElement::J1, "i", 160.0, 110.0),
                    node(BgElement::I, &format!("L={:.3}", p.inductance), 300.0, 30.0),
                    node(
                        BgElement::R,
                        &format!("R={:.3}", p.resistance),
                        300.0,
                        110.0,
                    ),
                    node(
                        BgElement::C,
                        &format!("C={:.4}", p.capacitance),
                        300.0,
                        190.0,
                    ),
                ];
                let bonds = vec![
                    BgBond { from: 0, to: 1 },
                    BgBond { from: 1, to: 2 },
                    BgBond { from: 1, to: 3 },
                    BgBond { from: 1, to: 4 },
                ];
                (nodes, bonds)
            }
            BondGraphPreset::DcMotor => {
                // Se:V — 1e — {R:Ra, I:La} — GY:Km — 1m — {I:J, R:b}
                let nodes = vec![
                    node(
                        BgElement::Se,
                        &format!("V={:.2}", p.motor_voltage),
                        20.0,
                        110.0,
                    ),
                    node(BgElement::J1, "i", 150.0, 110.0),
                    node(BgElement::R, &format!("Ra={:.3}", p.r_arm), 150.0, 30.0),
                    node(BgElement::I, &format!("La={:.3}", p.l_arm), 150.0, 190.0),
                    node(BgElement::GY, &format!("Km={:.3}", p.k_motor), 300.0, 110.0),
                    node(BgElement::J1, "\u{03C9}", 440.0, 110.0),
                    node(BgElement::I, &format!("J={:.3}", p.inertia), 440.0, 30.0),
                    node(BgElement::R, &format!("b={:.3}", p.friction), 440.0, 190.0),
                ];
                let bonds = vec![
                    BgBond { from: 0, to: 1 },
                    BgBond { from: 1, to: 2 },
                    BgBond { from: 1, to: 3 },
                    BgBond { from: 1, to: 4 },
                    BgBond { from: 4, to: 5 },
                    BgBond { from: 5, to: 6 },
                    BgBond { from: 5, to: 7 },
                ];
                (nodes, bonds)
            }
        }
    }

    /// Derive the linear **state-space** `dx/dt = A·x + B·u` for this preset
    /// from its bond-graph junction equations under **integral causality**. The
    /// state is the vector of **energy variables on the I and C elements**.
    ///
    /// **Mass-spring-damper** — state `x = [q, p]` (`q` = spring displacement on
    /// **C**, `p = m·v` = momentum on **I**); the 1-junction (common velocity
    /// `v = p/m`) gives `Σe = 0 → dp/dt = F − b·(p/m) − k·q` and `dq/dt = v`:
    /// ```text
    /// A = [[0,    1/m  ],   B = [0,   u = F
    ///      [-k,  -b/m ]]        1]
    /// ```
    /// Eliminating `p` reproduces the known ODE `m x'' + b x' + k x = F`.
    ///
    /// **Series RLC** — state `x = [q, λ]` (`q` = charge on **C**, `λ = L·i` =
    /// flux on **I**); the series 1-junction (common current `i = λ/L`) gives
    /// `dλ/dt = V − R·(λ/L) − q/C` and `dq/dt = i`:
    /// ```text
    /// A = [[0,    1/L  ],   B = [0,   u = V
    ///      [-1/C, -R/L]]        1]
    /// ```
    /// reproducing `L q'' + R q' + q/C = V`.
    ///
    /// **DC motor** — state `x = [λ, p]` (`λ = Lₐ·i` armature flux on **I**,
    /// `p = J·ω` rotor momentum on **I**); the GY couples the two 1-junctions
    /// (`e_back = Km·ω`, `τ = Km·i`):
    /// ```text
    /// dλ/dt = V − Rₐ·(λ/Lₐ) − Km·(p/J)
    /// dp/dt = Km·(λ/Lₐ) − b·(p/J)
    /// A = [[-Rₐ/Lₐ, -Km/J],   B = [1,   u = V
    ///      [ Km/Lₐ, -b/J ]]        0]
    /// ```
    /// the standard two-state DC-motor model.
    pub fn state_space(self, p: &BondGraphParams) -> StateSpace {
        match self {
            BondGraphPreset::MassSpringDamper => {
                let m = p.mass.max(1e-12);
                StateSpace {
                    a: vec![vec![0.0, 1.0 / m], vec![-p.stiffness, -p.damping / m]],
                    b: vec![0.0, 1.0],
                    u: p.force,
                    state_names: vec!["q (displacement)", "p (momentum)"],
                }
            }
            BondGraphPreset::Rlc => {
                let l = p.inductance.max(1e-12);
                let c = p.capacitance.max(1e-12);
                StateSpace {
                    a: vec![vec![0.0, 1.0 / l], vec![-1.0 / c, -p.resistance / l]],
                    b: vec![0.0, 1.0],
                    u: p.voltage,
                    state_names: vec!["q (charge)", "\u{03BB} (flux)"],
                }
            }
            BondGraphPreset::DcMotor => {
                let la = p.l_arm.max(1e-12);
                let j = p.inertia.max(1e-12);
                StateSpace {
                    a: vec![
                        vec![-p.r_arm / la, -p.k_motor / j],
                        vec![p.k_motor / la, -p.friction / j],
                    ],
                    b: vec![1.0, 0.0],
                    u: p.motor_voltage,
                    state_names: vec!["\u{03BB} (armature flux)", "p (rotor momentum)"],
                }
            }
        }
    }

    /// Undamped natural frequency `ωn` (rad/s) and damping ratio `ζ` of this
    /// preset's 2nd-order response, computed directly from the parameters (the
    /// closed-form the readout reports). Returns `None` for a degenerate (zero
    /// stiffness / capacitance) case where `ωn` is undefined.
    ///
    /// - MSD: `ωn = √(k/m)`, `ζ = b / (2√(k·m))`.
    /// - RLC: `ωn = 1/√(L·C)`, `ζ = (R/2)·√(C/L)`.
    /// - DC motor: characteristic poly `Lₐ J s² + (Rₐ J + b Lₐ) s + (Rₐ b + Km²)`
    ///   → `ωn = √((Rₐ b + Km²)/(Lₐ J))`, `ζ = (Rₐ J + b Lₐ)/(2√(Lₐ J (Rₐ b + Km²)))`.
    pub fn natural_freq_damping(self, p: &BondGraphParams) -> Option<(f64, f64)> {
        match self {
            BondGraphPreset::MassSpringDamper => {
                if p.stiffness <= 0.0 || p.mass <= 0.0 {
                    return None;
                }
                let wn = (p.stiffness / p.mass).sqrt();
                let zeta = p.damping / (2.0 * (p.stiffness * p.mass).sqrt());
                Some((wn, zeta))
            }
            BondGraphPreset::Rlc => {
                if p.capacitance <= 0.0 || p.inductance <= 0.0 {
                    return None;
                }
                let wn = 1.0 / (p.inductance * p.capacitance).sqrt();
                let zeta = 0.5 * p.resistance * (p.capacitance / p.inductance).sqrt();
                Some((wn, zeta))
            }
            BondGraphPreset::DcMotor => {
                let la = p.l_arm;
                let j = p.inertia;
                if la <= 0.0 || j <= 0.0 {
                    return None;
                }
                let k0 = p.r_arm * p.friction + p.k_motor * p.k_motor;
                if k0 <= 0.0 {
                    return None;
                }
                let wn = (k0 / (la * j)).sqrt();
                let zeta = (p.r_arm * j + p.friction * la) / (2.0 * (la * j * k0).sqrt());
                Some((wn, zeta))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// One solved run: the sampled `(t, x)` trajectory plus the state captions for
/// the plot legend.
#[derive(Clone, Debug, Default)]
pub struct BondGraphSolution {
    /// `(time, state-vector)` samples from the RK4 integration.
    pub samples: Vec<(f64, Vec<f64>)>,
    /// Caption per state variable.
    pub state_names: Vec<&'static str>,
}

/// Persistent state for the Bond Graph workbench: the chosen preset, its element
/// parameters and the latest solution.
pub struct BondGraphWorkbenchState {
    /// The active preset (drives the bond graph + the state-space).
    pub preset: BondGraphPreset,
    /// All element parameters (only the active preset's are read).
    pub params: BondGraphParams,
    /// The most recent solve result (empty until **Solve** runs).
    pub solution: BondGraphSolution,
}

impl Default for BondGraphWorkbenchState {
    fn default() -> Self {
        let mut s = Self {
            preset: BondGraphPreset::MassSpringDamper,
            params: BondGraphParams::default(),
            solution: BondGraphSolution::default(),
        };
        // Seed a solved trajectory so the panel shows a response immediately.
        s.solve();
        s
    }
}

impl BondGraphWorkbenchState {
    /// Derive the active preset's state-space and integrate it (RK4), storing
    /// the trajectory. Shared by the in-panel **Solve** button and the
    /// `bondgraph.solve` bridge id so both run the SAME path.
    pub fn solve(&mut self) {
        let ss = self.preset.state_space(&self.params);
        // ~600 samples across the window — smooth plot, cheap integration.
        let t_end = self.params.duration.max(1e-3);
        let dt = (t_end / 600.0).max(1e-5);
        let samples = ss.integrate_rk4(t_end, dt);
        self.solution = BondGraphSolution {
            samples,
            state_names: ss.state_names,
        };
    }

    /// Captions of every control the agent bridge can `SetControl`.
    pub fn agent_control_names() -> &'static [&'static str] {
        &[
            "Preset",
            "Mass",
            "Damping",
            "Stiffness",
            "Force",
            "Resistance",
            "Inductance",
            "Capacitance",
            "Voltage",
            "Armature resistance",
            "Armature inductance",
            "Rotor inertia",
            "Rotor friction",
            "Motor constant",
            "Motor voltage",
            "Duration",
        ]
    }

    /// Set one labelled control by caption for the agent `SetControl` bridge.
    /// Fail-loud on an unknown caption / wrong type / non-finite value; no state
    /// is written on error and nothing panics. **`Preset`** takes a preset id
    /// string (`msd` / `rlc` / `dcmotor`, plus aliases); every other control is
    /// a finite number written to the matching parameter.
    pub fn agent_set(
        &mut self,
        name: &str,
        value: &crate::agent_commands::AgentValue,
    ) -> Result<(), String> {
        // Helper: read a finite f64 or fail loud with the control name.
        let num = |value: &crate::agent_commands::AgentValue| -> Result<f64, String> {
            let v = value.as_f64()?;
            if !v.is_finite() {
                return Err(format!("{name}: value must be finite, got {v}"));
            }
            Ok(v)
        };
        match name {
            "Preset" => {
                let s = value.as_str()?;
                match BondGraphPreset::from_id(s) {
                    Some(pr) => {
                        self.preset = pr;
                        Ok(())
                    }
                    None => Err(format!(
                        "Preset: unknown id {s:?} (use msd / rlc / dcmotor)"
                    )),
                }
            }
            "Mass" => {
                self.params.mass = num(value)?;
                Ok(())
            }
            "Damping" => {
                self.params.damping = num(value)?;
                Ok(())
            }
            "Stiffness" => {
                self.params.stiffness = num(value)?;
                Ok(())
            }
            "Force" => {
                self.params.force = num(value)?;
                Ok(())
            }
            "Resistance" => {
                self.params.resistance = num(value)?;
                Ok(())
            }
            "Inductance" => {
                self.params.inductance = num(value)?;
                Ok(())
            }
            "Capacitance" => {
                self.params.capacitance = num(value)?;
                Ok(())
            }
            "Voltage" => {
                self.params.voltage = num(value)?;
                Ok(())
            }
            "Armature resistance" => {
                self.params.r_arm = num(value)?;
                Ok(())
            }
            "Armature inductance" => {
                self.params.l_arm = num(value)?;
                Ok(())
            }
            "Rotor inertia" => {
                self.params.inertia = num(value)?;
                Ok(())
            }
            "Rotor friction" => {
                self.params.friction = num(value)?;
                Ok(())
            }
            "Motor constant" => {
                self.params.k_motor = num(value)?;
                Ok(())
            }
            "Motor voltage" => {
                self.params.motor_voltage = num(value)?;
                Ok(())
            }
            "Duration" => {
                let v = num(value)?;
                if v <= 0.0 {
                    return Err(format!("Duration: must be > 0, got {v}"));
                }
                self.params.duration = v;
                Ok(())
            }
            other => Err(format!("unknown bondgraph control: {other:?}")),
        }
    }

    /// Readout for the agent `ReadReadout` bridge: the preset, the derived ODE
    /// order, the natural frequency & damping ratio, and the final state. Always
    /// `Some` (the workbench always has a preset).
    pub fn agent_readout(&self) -> Option<String> {
        let ss = self.preset.state_space(&self.params);
        let order = ss.order();
        let wd = self.preset.natural_freq_damping(&self.params);
        let wd_txt = match wd {
            Some((wn, zeta)) => {
                let regime = if zeta < 1.0 {
                    "underdamped"
                } else if (zeta - 1.0).abs() < 1e-9 {
                    "critically damped"
                } else {
                    "overdamped"
                };
                format!("\u{03C9}n={wn:.4} rad/s \u{00B7} \u{03B6}={zeta:.4} ({regime})")
            }
            None => "\u{03C9}n/\u{03B6} undefined (degenerate parameters)".to_string(),
        };
        let final_txt = match self.solution.samples.last() {
            Some((t, x)) => {
                let parts: Vec<String> = self
                    .solution
                    .state_names
                    .iter()
                    .zip(x.iter())
                    .map(|(n, v)| format!("{n}={v:.4}"))
                    .collect();
                format!("@t={t:.3}s: {}", parts.join(", "))
            }
            None => "(not solved)".to_string(),
        };
        Some(format!(
            "Bond graph \u{00B7} {} \u{00B7} ODE order {order} \u{00B7} {wd_txt} \u{00B7} {final_txt}",
            self.preset.label()
        ))
    }
}

// ---------------------------------------------------------------------------
// Bridge run action (solve)
// ---------------------------------------------------------------------------

/// Run the derive-then-integrate solve (the in-panel **Solve** action). Factored
/// out so the button and the `bondgraph.solve` bridge id share one path.
pub(crate) fn run(app: &mut ValenxApp) {
    app.bondgraph.solve();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the Bond Graph workbench. A no-op unless toggled on via View -> Bond
/// Graph.
pub fn draw_bondgraph_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_bondgraph_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_bondgraph_workbench",
        "Bond Graph (multi-domain systems modelling)",
        bondgraph_workbench_body,
    );
    if close {
        app.show_bondgraph_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn bondgraph_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "In-house bond-graph systems modeller [a bond graph models a physical system by power \
             flow: each bond carries an effort e and a flow f (e\u{00B7}f = power), and the same R / \
             C / I / Se / Sf / TF / GY / junction elements span mechanical, electrical, hydraulic \
             and thermal domains. Pick a preset, set its element parameters, and press Solve: the \
             standard bond-graph derivation builds dx/dt = A\u{00B7}x + B\u{00B7}u with the state = the \
             energy variables on the I and C elements, then it is integrated (RK4) and plotted. \
             Preset-based V1; general arbitrary-graph derivation later].",
        )
        .weak()
        .small(),
    );
    ui.separator();

    // --- Preset selector ----------------------------------------------------
    ui.horizontal(|ui| {
        let lbl = ui.label("Preset");
        egui::ComboBox::from_id_source("bondgraph_preset")
            .selected_text(app.bondgraph.preset.label())
            .show_ui(ui, |ui| {
                for pr in BondGraphPreset::ALL {
                    ui.selectable_value(&mut app.bondgraph.preset, pr, pr.label());
                }
            })
            .response
            .labelled_by(lbl.id);
        if ui
            .button("\u{25B6} Solve")
            .on_hover_text(
                "Derive the bond-graph state equations (dx/dt = A\u{00B7}x + B\u{00B7}u) and integrate \
                 them with RK4.",
            )
            .clicked()
        {
            app.bondgraph.solve();
        }
    });

    ui.separator();

    // --- Element parameter controls (per preset) ----------------------------
    ui.label(egui::RichText::new("Element parameters").strong());
    let preset = app.bondgraph.preset;
    let p = &mut app.bondgraph.params;
    egui::Grid::new("bondgraph_params")
        .num_columns(2)
        .spacing([12.0, 4.0])
        .show(ui, |ui| match preset {
            BondGraphPreset::MassSpringDamper => {
                param_row(ui, "Mass", "I: mass m (kg)", &mut p.mass, 0.01);
                param_row(
                    ui,
                    "Damping",
                    "R: damping b (N\u{00B7}s/m)",
                    &mut p.damping,
                    0.01,
                );
                param_row(
                    ui,
                    "Stiffness",
                    "C: stiffness k (N/m), compliance 1/k",
                    &mut p.stiffness,
                    0.1,
                );
                param_row(ui, "Force", "Se: applied force F (N)", &mut p.force, 0.1);
            }
            BondGraphPreset::Rlc => {
                param_row(
                    ui,
                    "Resistance",
                    "R: resistance (\u{03A9})",
                    &mut p.resistance,
                    0.01,
                );
                param_row(
                    ui,
                    "Inductance",
                    "I: inductance L (H)",
                    &mut p.inductance,
                    0.01,
                );
                param_row(
                    ui,
                    "Capacitance",
                    "C: capacitance (F)",
                    &mut p.capacitance,
                    0.0001,
                );
                param_row(
                    ui,
                    "Voltage",
                    "Se: source voltage V (V)",
                    &mut p.voltage,
                    0.1,
                );
            }
            BondGraphPreset::DcMotor => {
                param_row(
                    ui,
                    "Armature resistance",
                    "R: armature resistance Ra (\u{03A9})",
                    &mut p.r_arm,
                    0.01,
                );
                param_row(
                    ui,
                    "Armature inductance",
                    "I: armature inductance La (H)",
                    &mut p.l_arm,
                    0.01,
                );
                param_row(
                    ui,
                    "Rotor inertia",
                    "I: rotor inertia J (kg\u{00B7}m\u{00B2})",
                    &mut p.inertia,
                    0.001,
                );
                param_row(
                    ui,
                    "Rotor friction",
                    "R: rotor friction b (N\u{00B7}m\u{00B7}s)",
                    &mut p.friction,
                    0.001,
                );
                param_row(
                    ui,
                    "Motor constant",
                    "GY: motor / back-EMF constant Km (N\u{00B7}m/A)",
                    &mut p.k_motor,
                    0.001,
                );
                param_row(
                    ui,
                    "Motor voltage",
                    "Se: supply voltage V (V)",
                    &mut p.motor_voltage,
                    0.1,
                );
            }
        });
    ui.add_space(4.0);
    param_row_outside_grid(ui, "Duration", "simulated time (s)", &mut p.duration, 0.1);

    ui.separator();

    // --- Derived ODE summary ------------------------------------------------
    let ss = app.bondgraph.preset.state_space(&app.bondgraph.params);
    ui.label(
        egui::RichText::new(format!(
            "Derived state-space: ODE order {} (state = energy variables on the I/C elements)",
            ss.order()
        ))
        .strong(),
    );
    if let Some((wn, zeta)) = app
        .bondgraph
        .preset
        .natural_freq_damping(&app.bondgraph.params)
    {
        let regime = if zeta < 1.0 {
            "underdamped"
        } else if (zeta - 1.0).abs() < 1e-9 {
            "critically damped"
        } else {
            "overdamped"
        };
        ui.label(format!(
            "\u{03C9}n = {wn:.4} rad/s   ({:.4} Hz)   \u{00B7}   \u{03B6} = {zeta:.4}  ({regime})",
            wn / (2.0 * std::f64::consts::PI)
        ));
    }

    // --- The bond-graph canvas ----------------------------------------------
    ui.add_space(4.0);
    draw_bond_graph_canvas(app, ui);

    // --- Response plot ------------------------------------------------------
    ui.add_space(6.0);
    ui.label(egui::RichText::new("State response").strong());
    let sol = &app.bondgraph.solution;
    if sol.samples.len() < 2 {
        ui.label(egui::RichText::new("Press Solve to integrate the response.").weak());
    } else {
        let n_states = sol.state_names.len();
        managed_plot_mem_cfg(
            ui,
            "bondgraph_response",
            200.0,
            |plot| plot.legend(Legend::default()),
            |pui| {
                for s in 0..n_states {
                    let pts: PlotPoints = sol
                        .samples
                        .iter()
                        .map(|(t, x)| [*t, x.get(s).copied().unwrap_or(0.0)])
                        .collect();
                    pui.line(Line::new(pts).name(sol.state_names[s]));
                }
            },
        );
    }
}

/// One labelled `DragValue` parameter row inside the params grid. The caption is
/// a named label the DragValue is `labelled_by`, so the agent bridge / a screen
/// reader can find the spin button by its caption text (the AI-drivable name).
fn param_row(ui: &mut egui::Ui, caption: &str, hover: &str, value: &mut f64, speed: f64) {
    let lbl = ui.label(caption);
    ui.add(egui::DragValue::new(value).speed(speed).max_decimals(6))
        .labelled_by(lbl.id)
        .on_hover_text(hover);
    ui.end_row();
}

/// Same as [`param_row`] but for a row drawn outside the grid (horizontal).
fn param_row_outside_grid(
    ui: &mut egui::Ui,
    caption: &str,
    hover: &str,
    value: &mut f64,
    speed: f64,
) {
    ui.horizontal(|ui| {
        let lbl = ui.label(caption);
        ui.add(egui::DragValue::new(value).speed(speed).max_decimals(6))
            .labelled_by(lbl.id)
            .on_hover_text(hover);
    });
}

// ---------------------------------------------------------------------------
// Bond-graph canvas (reuses the in-house node-graph drawing idiom)
// ---------------------------------------------------------------------------

/// Fixed node geometry (canvas-local points).
const BG_NODE_W: f32 = 96.0;
const BG_NODE_H: f32 = 44.0;

/// Draw the active preset's bond graph: element nodes joined by power bonds.
/// Read-only (no interaction in V1 — the graph is the preset's; the user edits
/// the system through the parameter controls), so this is pure painting using
/// the same bezier-bond / node-rect idiom as the node-graph workbench.
fn draw_bond_graph_canvas(app: &mut ValenxApp, ui: &mut egui::Ui) {
    let (nodes, bonds) = app.bondgraph.preset.bond_graph(&app.bondgraph.params);

    let desired = egui::vec2(ui.available_width(), 240.0);
    let (resp, painter) = ui.allocate_painter(desired, egui::Sense::hover());
    let origin = resp.rect.min.to_vec2();

    painter.rect_filled(resp.rect, 4.0, egui::Color32::from_rgb(24, 26, 32));
    painter.rect_stroke(
        resp.rect,
        4.0,
        egui::Stroke::new(1.0, egui::Color32::from_rgb(50, 55, 68)),
    );

    // Centre of a node (for bond endpoints).
    let center = |n: &BgNode| -> egui::Pos2 {
        egui::pos2(
            n.pos.x + origin.x + BG_NODE_W * 0.5,
            n.pos.y + origin.y + BG_NODE_H * 0.5,
        )
    };

    // --- Power bonds (drawn first, under the nodes) -------------------------
    for b in &bonds {
        if let (Some(a), Some(c)) = (nodes.get(b.from), nodes.get(b.to)) {
            let pa = center(a);
            let pc = center(c);
            painter.line_segment(
                [pa, pc],
                egui::Stroke::new(2.0, egui::Color32::from_rgb(120, 170, 255)),
            );
            // A small half-arrow at the destination marks the bond's power
            // sign-convention end (the bond-graph half-arrow), drawn as a short
            // tick toward the source.
            let dir = (pa - pc).normalized();
            let tick = pc + dir * 12.0;
            let perp = egui::vec2(-dir.y, dir.x) * 4.0;
            painter.line_segment(
                [pc, tick + perp],
                egui::Stroke::new(2.0, egui::Color32::from_rgb(120, 170, 255)),
            );
        }
    }

    // --- Element nodes ------------------------------------------------------
    for n in &nodes {
        let tl = egui::pos2(n.pos.x + origin.x, n.pos.y + origin.y);
        let body = egui::Rect::from_min_size(tl, egui::vec2(BG_NODE_W, BG_NODE_H));
        painter.rect_filled(body, 5.0, egui::Color32::from_rgb(40, 44, 54));
        painter.rect_stroke(body, 5.0, egui::Stroke::new(1.5, n.element.color()));
        // Element symbol (top) + caption (bottom).
        painter.text(
            egui::pos2(body.center().x, tl.y + 13.0),
            egui::Align2::CENTER_CENTER,
            n.element.symbol(),
            egui::FontId::proportional(14.0),
            egui::Color32::from_gray(235),
        );
        painter.text(
            egui::pos2(body.center().x, tl.y + BG_NODE_H - 12.0),
            egui::Align2::CENTER_CENTER,
            &n.caption,
            egui::FontId::monospace(10.0),
            egui::Color32::from_gray(185),
        );
    }
}

// ---------------------------------------------------------------------------
// Tests (unit)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::agent_commands::AgentValue;

    /// Fit a 2nd-order step response `x(t)` and confirm it satisfies the target
    /// ODE residual to tolerance by finite differences. Returns the RMS residual
    /// of `a·x'' + b·x' + c·x − rhs` over the interior samples.
    fn ode_residual_rms(samples: &[(f64, f64)], a: f64, b: f64, c: f64, rhs: f64) -> f64 {
        let mut acc = 0.0;
        let mut cnt = 0usize;
        for i in 1..samples.len() - 1 {
            let (t0, x0) = samples[i - 1];
            let (_t1, x1) = samples[i];
            let (t2, x2) = samples[i + 1];
            let dt = (t2 - t0) * 0.5;
            let xp = (x2 - x0) / (2.0 * dt);
            let xpp = (x2 - 2.0 * x1 + x0) / (dt * dt);
            let r = a * xpp + b * xp + c * x0 - rhs;
            acc += r * r;
            cnt += 1;
            let _ = t0;
        }
        if cnt == 0 {
            return 0.0;
        }
        (acc / cnt as f64).sqrt()
    }

    #[test]
    fn msd_state_space_matches_analytic_ode() {
        // CANONICAL #1: the mass-spring-damper state-space must reproduce
        // m x'' + b x' + k x = F. We integrate dx/dt = A x + B u (state [q,p])
        // and check the *displacement* q against the analytic ODE residual.
        // A well-damped case (zeta = 8/(2*sqrt(15*2)) ~ 0.73) so the transient
        // settles well within the window and the steady-state check is exact.
        let p = BondGraphParams {
            mass: 2.0,
            damping: 8.0,
            stiffness: 15.0,
            force: 3.0,
            ..Default::default()
        };
        let ss = BondGraphPreset::MassSpringDamper.state_space(&p);
        assert_eq!(ss.order(), 2, "two energy stores -> 2nd-order ODE");
        let traj = ss.integrate_rk4(12.0, 1e-3);
        // Extract q (state index 0).
        let q: Vec<(f64, f64)> = traj.iter().map(|(t, x)| (*t, x[0])).collect();
        // The ODE residual (independent of damping) is the real validation that
        // the derived state-space IS m x'' + b x' + k x = F.
        let res = ode_residual_rms(&q, p.mass, p.damping, p.stiffness, p.force);
        assert!(
            res < 1e-2,
            "m x'' + b x' + k x = F residual too large: {res}"
        );
        // Steady state q_inf = F/k (the transient has fully decayed by 12 s).
        let q_inf = q.last().unwrap().1;
        assert!(
            (q_inf - p.force / p.stiffness).abs() < 1e-4,
            "steady displacement should approach F/k = {}, got {q_inf}",
            p.force / p.stiffness
        );
    }

    #[test]
    fn rlc_state_space_matches_analytic_ode() {
        // CANONICAL #2: the series-RLC state-space must reproduce
        // L q'' + R q' + q/C = V. State [q, lambda]; check charge q.
        let p = BondGraphParams {
            resistance: 3.0,
            inductance: 0.5,
            capacitance: 0.02,
            voltage: 5.0,
            ..Default::default()
        };
        let ss = BondGraphPreset::Rlc.state_space(&p);
        assert_eq!(ss.order(), 2);
        let traj = ss.integrate_rk4(2.0, 5e-4);
        let q: Vec<(f64, f64)> = traj.iter().map(|(t, x)| (*t, x[0])).collect();
        let res = ode_residual_rms(
            &q,
            p.inductance,
            p.resistance,
            1.0 / p.capacitance,
            p.voltage,
        );
        assert!(
            res < 1e-2,
            "L q'' + R q' + q/C = V residual too large: {res}"
        );
        // Steady charge q_inf = C*V (capacitor charges to source voltage).
        let q_inf = q.last().unwrap().1;
        assert!(
            (q_inf - p.capacitance * p.voltage).abs() < 1e-3,
            "steady charge should approach C*V = {}, got {q_inf}",
            p.capacitance * p.voltage
        );
    }

    #[test]
    fn msd_natural_freq_and_damping_are_correct() {
        // ωn = sqrt(k/m), ζ = b/(2 sqrt(k m)).
        let p = BondGraphParams {
            mass: 4.0,
            damping: 2.0,
            stiffness: 100.0,
            ..Default::default()
        };
        let (wn, zeta) = BondGraphPreset::MassSpringDamper
            .natural_freq_damping(&p)
            .unwrap();
        assert!((wn - 5.0).abs() < 1e-9, "wn = sqrt(100/4) = 5, got {wn}");
        assert!(
            (zeta - 2.0 / (2.0 * (100.0_f64 * 4.0).sqrt())).abs() < 1e-9,
            "zeta mismatch: {zeta}"
        );
    }

    #[test]
    fn rlc_natural_freq_and_damping_are_correct() {
        // ωn = 1/sqrt(LC), ζ = (R/2) sqrt(C/L).
        let p = BondGraphParams {
            resistance: 4.0,
            inductance: 1.0,
            capacitance: 0.0025,
            ..Default::default()
        };
        let (wn, zeta) = BondGraphPreset::Rlc.natural_freq_damping(&p).unwrap();
        assert!(
            (wn - 20.0).abs() < 1e-9,
            "wn = 1/sqrt(1*0.0025) = 20, got {wn}"
        );
        assert!(
            (zeta - 2.0 * (0.0025_f64 / 1.0).sqrt()).abs() < 1e-9,
            "zeta mismatch: {zeta}"
        );
    }

    #[test]
    fn dc_motor_reaches_known_steady_state_speed() {
        // The DC-motor steady-state angular velocity from
        // dλ/dt = dp/dt = 0:  ω∞ = Km·V / (Ra·b + Km²).
        let p = BondGraphParams {
            r_arm: 1.0,
            l_arm: 0.5,
            inertia: 0.01,
            friction: 0.1,
            k_motor: 0.1,
            motor_voltage: 12.0,
            duration: 6.0,
            ..Default::default()
        };
        let ss = BondGraphPreset::DcMotor.state_space(&p);
        assert_eq!(ss.order(), 2);
        let traj = ss.integrate_rk4(6.0, 1e-3);
        // omega = p / J (state index 1 is p).
        let (_t, x) = traj.last().unwrap();
        let omega = x[1] / p.inertia;
        let expected = p.k_motor * p.motor_voltage / (p.r_arm * p.friction + p.k_motor * p.k_motor);
        assert!(
            (omega - expected).abs() < 1e-2,
            "DC-motor steady speed should be Km*V/(Ra*b+Km^2) = {expected}, got {omega}"
        );
    }

    #[test]
    fn rk4_returns_initial_sample_for_bad_dt() {
        let ss = BondGraphPreset::MassSpringDamper.state_space(&BondGraphParams::default());
        assert_eq!(
            ss.integrate_rk4(1.0, 0.0).len(),
            1,
            "non-positive dt -> just x0"
        );
        assert_eq!(
            ss.integrate_rk4(-1.0, 0.1).len(),
            1,
            "non-positive t_end -> just x0"
        );
        assert_eq!(
            ss.integrate_rk4(f64::NAN, 0.1).len(),
            1,
            "non-finite t_end -> just x0"
        );
    }

    #[test]
    fn bond_graph_has_expected_elements_per_preset() {
        let p = BondGraphParams::default();
        // MSD: Se, 1, I, R, C (5 nodes, 4 bonds).
        let (nodes, bonds) = BondGraphPreset::MassSpringDamper.bond_graph(&p);
        assert_eq!(nodes.len(), 5);
        assert_eq!(bonds.len(), 4);
        assert!(nodes.iter().any(|n| n.element == BgElement::I));
        assert!(nodes.iter().any(|n| n.element == BgElement::C));
        assert!(nodes.iter().any(|n| n.element == BgElement::J1));
        // DC motor uses a gyrator.
        let (mnodes, _) = BondGraphPreset::DcMotor.bond_graph(&p);
        assert!(mnodes.iter().any(|n| n.element == BgElement::GY));
        // Two I elements (La and J).
        assert_eq!(
            mnodes.iter().filter(|n| n.element == BgElement::I).count(),
            2
        );
    }

    #[test]
    fn solve_populates_a_trajectory() {
        let mut s = BondGraphWorkbenchState::default();
        assert!(
            s.solution.samples.len() > 100,
            "default state seeds a solved trajectory"
        );
        assert_eq!(s.solution.state_names.len(), 2);
        // Re-solving after a parameter change keeps a trajectory.
        s.params.stiffness = 50.0;
        s.solve();
        assert!(s.solution.samples.len() > 100);
    }

    #[test]
    fn agent_set_preset_switches_system() {
        let mut s = BondGraphWorkbenchState::default();
        s.agent_set("Preset", &AgentValue::Str("rlc".to_string()))
            .expect("switch to RLC");
        assert_eq!(s.preset, BondGraphPreset::Rlc);
        s.agent_set("Preset", &AgentValue::Str("dcmotor".to_string()))
            .expect("switch to DC motor");
        assert_eq!(s.preset, BondGraphPreset::DcMotor);
        // Unknown preset is fail-loud.
        assert!(s
            .agent_set("Preset", &AgentValue::Str("frobnicate".to_string()))
            .is_err());
    }

    #[test]
    fn agent_set_numeric_params_round_trip() {
        let mut s = BondGraphWorkbenchState::default();
        s.agent_set("Mass", &AgentValue::Float(3.5)).unwrap();
        assert_eq!(s.params.mass, 3.5);
        s.agent_set("Stiffness", &AgentValue::Float(42.0)).unwrap();
        assert_eq!(s.params.stiffness, 42.0);
        s.agent_set("Capacitance", &AgentValue::Float(0.001))
            .unwrap();
        assert_eq!(s.params.capacitance, 0.001);
        // Non-finite rejected.
        assert!(s.agent_set("Mass", &AgentValue::Float(f64::NAN)).is_err());
        // Non-positive duration rejected.
        assert!(s.agent_set("Duration", &AgentValue::Float(0.0)).is_err());
        // Unknown control rejected.
        assert!(s.agent_set("Nope", &AgentValue::Float(1.0)).is_err());
    }

    #[test]
    fn control_names_listed_and_nonempty() {
        let names = BondGraphWorkbenchState::agent_control_names();
        assert!(names.contains(&"Preset"));
        assert!(names.contains(&"Mass"));
        assert!(names.contains(&"Duration"));
        assert!(names.contains(&"Motor constant"));
    }

    #[test]
    fn readout_reports_order_freq_and_final_state() {
        let s = BondGraphWorkbenchState::default();
        let r = s.agent_readout().expect("readout always present");
        assert!(r.contains("ODE order 2"), "got: {r}");
        assert!(r.contains("Mass-spring-damper"), "got: {r}");
        assert!(
            r.contains("\u{03C9}n="),
            "natural frequency in readout; got: {r}"
        );
        assert!(
            r.contains("\u{03B6}="),
            "damping ratio in readout; got: {r}"
        );
    }

    #[test]
    fn run_bridge_helper_solves_through_app() {
        let mut app = ValenxApp::default();
        app.bondgraph
            .agent_set("Preset", &AgentValue::Str("rlc".to_string()))
            .unwrap();
        app.bondgraph
            .agent_set("Duration", &AgentValue::Float(1.0))
            .unwrap();
        run(&mut app);
        assert!(app.bondgraph.solution.samples.len() > 100);
        // The readout reflects the RLC preset after the bridge solve.
        let r = app.bondgraph.agent_readout().unwrap();
        assert!(r.contains("Series RLC"), "got: {r}");
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
            draw_bondgraph_workbench(app, ctx);
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
        assert!(!app.show_bondgraph_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_bondgraph_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown() {
        let mut app = ValenxApp::default();
        app.show_bondgraph_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);
        assert!(!nodes.is_empty(), "a shown workbench produces a11y nodes");
    }

    #[test]
    fn numeric_controls_are_named_and_associated() {
        // Every parameter DragValue is a SpinButton and must be `labelled_by`
        // its caption so an AI / screen reader can find it by caption text.
        let mut app = ValenxApp::default();
        app.show_bondgraph_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            !spin_buttons.is_empty(),
            "expected the parameter numeric controls as spin buttons"
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
        // The default preset is the mass-spring-damper, so its captions appear.
        assert!(
            has_named_node(&nodes, "Mass"),
            "'Mass' caption is a named node"
        );
        assert!(
            has_named_node(&nodes, "Stiffness"),
            "'Stiffness' caption is a named node"
        );
        assert!(
            has_named_node(&nodes, "Duration"),
            "'Duration' caption is a named node"
        );
    }
}
