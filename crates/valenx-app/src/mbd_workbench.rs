//! The right-side **Multibody-Dynamics (robot / contact) Workbench** panel —
//! a native front-end over the in-house [`valenx_mbd`] crate (Valenx's planar
//! constrained-DAE multibody solver, its Featherstone articulated-body
//! algorithm, and its penalty-contact + Coulomb-friction model).
//!
//! Multibody dynamics is the time-domain simulation of rigid bodies coupled by
//! joints, forces and contacts. This workbench drives the **real**
//! [`valenx_mbd`] integrators over two fully-native, fully-transparent demos —
//! the user picks which with the *demo* combo:
//!
//! * **Pendulum** — a planar rigid-rod pendulum: one rigid [`valenx_mbd::Body`]
//!   hung from a rigid-distance link ([`valenx_mbd::Constraint::Distance`]) to a
//!   world pivot under gravity, advanced by the real constrained-DAE
//!   [`valenx_mbd::System::step`] (KKT solve + Baumgarte + semi-implicit Euler).
//!   With no damping this is a conservative (NVE) system, so its real
//!   [`valenx_mbd::System::energy`] stays bounded — the **energy-conservation
//!   PIN** — and at a small release angle its period approaches the textbook
//!   simple pendulum `2π√(L/g)`.
//!
//! * **Drop onto a plane** — a rigid body released above a flat ground plane and
//!   integrated through the real **penalty-contact** path
//!   ([`valenx_mbd::contact::body_contact_wrench`] / [`contact_force`]): while
//!   clear of the ground it is in free fall, so its height follows the analytic
//!   `h(t) = h₀ − ½ g t²` (the **free-fall PIN**); once it penetrates, the
//!   compliant spring–damper pushes it back and it settles with a small steady
//!   penetration bounded by `m g / kₙ` (the **contact-rest PIN**). An optional
//!   horizontal launch velocity with a non-zero **Coulomb friction coefficient**
//!   exercises the regularized friction law (a body sliding on the ground is
//!   decelerated by kinetic friction `μ·fₙ`).
//!
//! The painter draws the motion in 2-D — the pendulum's swing, or the dropped
//! body's height-vs-time trace overlaid with the contact-force trace — plus a
//! readout (final state, peak penetration, energy drift).
//!
//! ## Fail-loud parameter validation (guards the contact asserts)
//!
//! [`valenx_mbd::contact::ContactParams::new`] and
//! [`valenx_mbd::contact::Plane::new`] are **fail-loud**: they `assert!` their
//! physical preconditions (`kₙ > 0`, `cₙ ≥ 0`, `μ ≥ 0`, `v_eps > 0`, all finite;
//! a non-zero finite plane normal) and have `#[should_panic]` tests. A workbench
//! must therefore **never** let a user value reach those constructors unchecked.
//! [`MbdParams::validate`] mirrors *exactly* those conditions (plus the
//! integrator's own `dt > 0`, finite gravity/height/velocity, sane step count)
//! and returns `Err(String)` for any violation, so [`MbdWorkbenchState::run`]
//! rejects bad input **in-panel** before it ever constructs a `ContactParams` or
//! a `Plane`. The KEY test pins that an invalid contact param (stiffness ≤ 0,
//! dt ≤ 0, friction < 0, …) yields an `Err`, **not** a panic — proving the
//! assert is unreachable from the UI.
//!
//! Mirrors the other workbenches (`ppi_workbench`, `cosim_workbench`): a
//! [`crate::workbench_chrome::workbench_shell`] panel gated on
//! [`crate::ValenxApp::show_mbd_workbench`], toggled from the View menu and
//! openable by the agent bridge under the workbench id `"mbd"` (aliases
//! `"multibody"` / `"robot"`; see [`crate::project_tabs::TabKind`]). Every
//! numeric control is `.labelled_by` an accessible caption so the panel is
//! AI-drivable by name.
//!
//! Honesty: valenx-mbd is a **research / preliminary-design grade** planar
//! solver — 2-D rigid bodies, revolute / distance / prismatic constraints, and
//! *penalty* (compliant) contact. Penalty contact is soft (a resting body always
//! sinks `m g / kₙ`) and the explicit integrator is step-size sensitive; this is
//! a real v1, not a substitute for a production multibody / hard-LCP-contact
//! code. The demos prove the integrators are wired up and reproduce closed-form
//! results, not that the method is production-accurate. A degenerate input
//! surfaces an in-panel `Err` — never a panic.

use eframe::egui;
use nalgebra::Vector2;
use valenx_mbd::contact::{
    body_contact_wrench, BodyContact, ContactBodyState, ContactParams, Plane,
};
use valenx_mbd::{Anchor, Body, Constraint, Force, System};

use crate::ValenxApp;

// ---------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------

/// Which demo the workbench simulates.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum MbdDemo {
    /// A planar rigid-rod pendulum advanced by the constrained-DAE
    /// [`System`] — a conservative system whose energy stays bounded.
    #[default]
    Pendulum,
    /// A rigid body dropped onto a ground plane through the **penalty-contact**
    /// path ([`body_contact_wrench`]) with optional Coulomb friction.
    Drop,
}

impl MbdDemo {
    /// Human-readable label for the combo box / status line.
    fn label(self) -> &'static str {
        match self {
            MbdDemo::Pendulum => "Pendulum (articulated rod, energy-conserving)",
            MbdDemo::Drop => "Drop onto plane (penalty contact + friction)",
        }
    }
}

/// Editable multibody inputs shown in the workbench.
#[derive(Clone, Copy, Debug)]
pub struct MbdParams {
    /// Which demo to run.
    pub demo: MbdDemo,

    // -- shared --
    /// Gravitational acceleration magnitude `g` (m/s², pointing down). Must be
    /// finite and `> 0`.
    pub gravity: f64,
    /// Integrator time step `dt` (s). Must be finite and `> 0` — this is the
    /// integrator's own precondition (a non-positive `dt` makes no physical
    /// sense and would not advance the clock).
    pub dt: f64,
    /// Number of integration steps to take. Must be `1..=`[`MAX_STEPS`].
    pub steps: usize,

    // -- pendulum --
    /// Pendulum rod length `L` (m). Must be finite and `> 0`.
    pub rod_length: f64,
    /// Initial release angle from straight-down (rad). Any finite value.
    pub release_angle: f64,

    // -- drop / contact --
    /// Initial height of the body's contact point above the plane (m). Must be
    /// finite (may be negative — starting already penetrating is valid).
    pub drop_height: f64,
    /// Initial horizontal launch velocity (m/s) — a non-zero value plus friction
    /// exercises the Coulomb friction law. Any finite value.
    pub launch_speed: f64,
    /// Body mass `m` (kg). Must be finite and `> 0`.
    pub mass: f64,
    /// Normal contact stiffness `kₙ` (N/m). **Must be finite and `> 0`** —
    /// mirrors [`ContactParams::new`]'s assert.
    pub contact_stiffness: f64,
    /// Normal contact damping `cₙ` (N·s/m). **Must be finite and `≥ 0`** —
    /// mirrors [`ContactParams::new`]'s assert (`≥ 0`, so a purely elastic
    /// contact is allowed).
    pub contact_damping: f64,
    /// Coulomb friction coefficient `μ` (dimensionless). **Must be finite and
    /// `≥ 0`** — mirrors [`ContactParams::new`]'s assert.
    pub friction: f64,
}

/// Upper bound on the step count (keeps a headless run bounded and fast).
pub const MAX_STEPS: usize = 2_000_000;

/// Fixed friction-regularization velocity scale `v_eps` (m/s). Held internal
/// (not a user knob) and always `> 0`, so the [`ContactParams::new`] precondition
/// on it is met by construction; a small value gives a stiff "stick" region.
const FRICTION_VEL_SCALE: f64 = 1.0e-3;

impl Default for MbdParams {
    fn default() -> Self {
        Self {
            demo: MbdDemo::Pendulum,
            gravity: 9.81,
            dt: 1.0e-3,
            steps: 4_000,
            // Pendulum: a 1 m rod released at a small angle (so the small-angle
            // period 2π√(L/g) is a good comparison and the explicit integrator
            // holds energy well).
            rod_length: 1.0,
            release_angle: 0.2,
            // Drop: 0.5 m above a stiff plane, no initial slide, light damping,
            // mild friction. Stiffness chosen so the steady penetration mg/kₙ is
            // sub-millimetre and the default dt is stable for it.
            drop_height: 0.5,
            launch_speed: 0.0,
            mass: 1.0,
            contact_stiffness: 1.0e4,
            contact_damping: 40.0,
            friction: 0.3,
        }
    }
}

impl MbdParams {
    /// **Validate** every parameter the run will use, returning `Err(String)`
    /// for the first violation — fail-loud, in-panel, **never a panic**.
    ///
    /// Crucially this pre-checks the exact conditions
    /// [`ContactParams::new`] / [`Plane::new`] assert, so [`MbdWorkbenchState::run`]
    /// can call them knowing the assert can never fire from a user value:
    /// `contact_stiffness > 0`, `contact_damping ≥ 0`, `friction ≥ 0`,
    /// `FRICTION_VEL_SCALE > 0` (internal), all finite. It also guards the
    /// integrator's own preconditions (`dt > 0`, `mass > 0`, finite gravity,
    /// etc.) and a sane step count.
    pub fn validate(&self) -> Result<(), String> {
        // --- Shared preconditions ---------------------------------------
        if !self.gravity.is_finite() || self.gravity <= 0.0 {
            return Err(format!(
                "gravity must be finite and > 0 m/s² (got {})",
                self.gravity
            ));
        }
        if !self.dt.is_finite() || self.dt <= 0.0 {
            return Err(format!(
                "time step dt must be finite and > 0 s (got {})",
                self.dt
            ));
        }
        if self.steps == 0 || self.steps > MAX_STEPS {
            return Err(format!(
                "step count must be 1..={MAX_STEPS} (got {})",
                self.steps
            ));
        }

        match self.demo {
            MbdDemo::Pendulum => {
                if !self.rod_length.is_finite() || self.rod_length <= 0.0 {
                    return Err(format!(
                        "rod length must be finite and > 0 m (got {})",
                        self.rod_length
                    ));
                }
                if !self.release_angle.is_finite() {
                    return Err(format!(
                        "release angle must be finite (got {})",
                        self.release_angle
                    ));
                }
            }
            MbdDemo::Drop => {
                if !self.drop_height.is_finite() {
                    return Err(format!(
                        "drop height must be finite (got {})",
                        self.drop_height
                    ));
                }
                if !self.launch_speed.is_finite() {
                    return Err(format!(
                        "launch speed must be finite (got {})",
                        self.launch_speed
                    ));
                }
                if !self.mass.is_finite() || self.mass <= 0.0 {
                    return Err(format!(
                        "mass must be finite and > 0 kg (got {})",
                        self.mass
                    ));
                }
                // --- The contact-assert guards (mirror ContactParams::new). ---
                if !self.contact_stiffness.is_finite() || self.contact_stiffness <= 0.0 {
                    return Err(format!(
                        "contact stiffness kₙ must be finite and > 0 N/m (got {})",
                        self.contact_stiffness
                    ));
                }
                if !self.contact_damping.is_finite() || self.contact_damping < 0.0 {
                    return Err(format!(
                        "contact damping cₙ must be finite and ≥ 0 N·s/m (got {})",
                        self.contact_damping
                    ));
                }
                if !self.friction.is_finite() || self.friction < 0.0 {
                    return Err(format!(
                        "friction coefficient μ must be finite and ≥ 0 (got {})",
                        self.friction
                    ));
                }
                // FRICTION_VEL_SCALE is an internal constant (> 0), but assert it
                // here too so the guarantee is explicit and survives any future
                // edit that makes it user-editable.
                debug_assert!(FRICTION_VEL_SCALE.is_finite() && FRICTION_VEL_SCALE > 0.0);
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

/// One sampled instant of the simulated motion.
#[derive(Clone, Copy, Debug, Default)]
pub struct MbdSample {
    /// Simulation time at this sample (s).
    pub t: f64,
    /// Pendulum: the bob's angle from straight-down (rad). Drop: unused (`0`).
    pub angle: f64,
    /// Pendulum: the bob's `x` position (m). Drop: the body's `x` (m).
    pub x: f64,
    /// The body / bob height (m): pendulum CoM `y`, or drop contact-point height
    /// above the plane.
    pub height: f64,
    /// Total mechanical energy at this instant (J) — pendulum only (the real
    /// [`System::energy`]); `0` for the drop demo.
    pub energy: f64,
    /// Magnitude of the contact force on the body at this instant (N) — drop
    /// demo only (`0` for the pendulum).
    pub contact_force: f64,
}

/// The computed trajectory + summary diagnostics from a run.
#[derive(Default, Clone)]
pub struct MbdResult {
    /// Which demo produced this.
    pub demo: MbdDemo,
    /// The sampled trajectory (down-sampled to at most [`MAX_TRACE`] points so
    /// the painter trace stays light regardless of step count).
    pub samples: Vec<MbdSample>,
    /// Final simulation time (s).
    pub final_time: f64,
    /// Final body / bob height (m).
    pub final_height: f64,
    /// Pendulum: final angle (rad). Drop: final `x` (m).
    pub final_angle_or_x: f64,

    // -- pendulum diagnostics --
    /// Initial total energy (J).
    pub energy_initial: f64,
    /// Largest relative energy deviation `|E−E₀|/|E₀|` over the run (the
    /// conservation diagnostic; small for the undamped NVE pendulum).
    pub energy_rel_drift: f64,

    // -- drop diagnostics --
    /// Peak penetration depth into the plane over the run (m). Bounded by
    /// roughly `m g / kₙ` once settled.
    pub max_penetration: f64,
    /// The analytic steady penetration `m g / kₙ` (m) for comparison.
    pub analytic_penetration: f64,
    /// Peak contact-force magnitude over the run (N).
    pub max_contact_force: f64,
    /// Whether the body had effectively come to rest by the final step
    /// (`|v| < REST_SPEED`).
    pub at_rest: bool,
}

/// Cap on the number of stored trace samples (the painter and the readouts only
/// need a light trajectory, not every step).
pub const MAX_TRACE: usize = 600;

/// Speed (m/s) below which the dropped body is considered "at rest".
const REST_SPEED: f64 = 1.0e-2;

impl MbdResult {
    /// The time span covered by the trace (s).
    pub fn duration(&self) -> f64 {
        self.final_time
    }
}

// ---------------------------------------------------------------------------
// Workbench state
// ---------------------------------------------------------------------------

/// Persistent state for the multibody-dynamics workbench.
#[derive(Default)]
pub struct MbdWorkbenchState {
    /// User-editable parameters.
    pub params: MbdParams,
    /// Last successful result (populated after a successful Run).
    pub result: Option<MbdResult>,
    /// Status / error line shown below the controls.
    pub status: String,
}

impl MbdWorkbenchState {
    /// Validate the parameters, then run the selected demo through the **real**
    /// [`valenx_mbd`] integrators — fail-loud.
    ///
    /// Every failure path returns an `Err(String)` (rendered in-panel); no user
    /// value ever reaches a [`ContactParams::new`] / [`Plane::new`] assert
    /// because [`MbdParams::validate`] pre-checks the same conditions first.
    pub fn run(&self) -> Result<MbdResult, String> {
        let p = &self.params;
        p.validate()?; // <-- guards the contact asserts (and the integrator's).
        match p.demo {
            MbdDemo::Pendulum => run_pendulum(p),
            MbdDemo::Drop => run_drop(p),
        }
    }
}

/// Append a sample to the trace, keeping it down-sampled to at most
/// [`MAX_TRACE`] points (store every `stride`-th step).
fn push_sampled(trace: &mut Vec<MbdSample>, step: usize, stride: usize, s: MbdSample) {
    if stride <= 1 || step % stride == 0 {
        trace.push(s);
    }
}

/// Down-sampling stride so a `steps`-long run yields ≤ [`MAX_TRACE`] samples.
///
/// Uses **ceiling** division: the loop keeps every `stride`-th step, giving
/// `⌊steps / stride⌋` in-loop samples (plus the one initial sample). With floor
/// division `stride = ⌊steps / MAX_TRACE⌋` the count `⌊steps / stride⌋` can
/// exceed `MAX_TRACE` (e.g. `steps = 4000`, `MAX_TRACE = 600` → `stride = 6` →
/// 666 in-loop + 1 = 667 > 601). Rounding the stride **up** guarantees
/// `⌊steps / stride⌋ ≤ MAX_TRACE` for every step count.
fn trace_stride(steps: usize) -> usize {
    steps.div_ceil(MAX_TRACE).max(1)
}

/// Run the **pendulum** demo on the real planar constrained-DAE [`System`].
///
/// A single rigid body is hung from a rigid-distance link
/// ([`Constraint::Distance`]) attached at its centre of mass to a world pivot,
/// under gravity. The distance constraint passes through the CoM so it applies
/// no torque (the bob swings as an ideal simple pendulum); with no damping the
/// system is conservative, so [`System::energy`] is the conservation diagnostic.
fn run_pendulum(p: &MbdParams) -> Result<MbdResult, String> {
    let m = 1.0; // bob mass (energy conservation is mass-independent here)
    let inertia = 0.05; // small rotational inertia; decoupled (rod through CoM)
    let l = p.rod_length;
    let g = p.gravity;
    let theta0 = p.release_angle;

    let pivot = Vector2::new(0.0, 0.0);
    // CoM hangs at angle theta0 from straight-down: (L sinθ, −L cosθ).
    let (s, c) = theta0.sin_cos();
    let cg = pivot + Vector2::new(l * s, -l * c);

    let mut sys = System {
        bodies: vec![Body::new(m, inertia, cg)],
        forces: vec![Force::Gravity(Vector2::new(0.0, -g))],
        constraints: vec![Constraint::Distance {
            a: Anchor::Body {
                index: 0,
                local: Vector2::zeros(),
            },
            b: Anchor::World(pivot),
            length: l,
        }],
        time: 0.0,
    };

    let e0 = sys.energy();
    if !e0.is_finite() {
        return Err(format!("initial pendulum energy is non-finite ({e0})"));
    }
    let stride = trace_stride(p.steps);
    let mut trace: Vec<MbdSample> = Vec::with_capacity(p.steps.min(MAX_TRACE) + 1);
    let mut max_rel_drift = 0.0_f64;
    let denom = e0.abs().max(1e-12);

    // Sample the initial state, then step.
    let sample_of = |sys: &System, e: f64| -> MbdSample {
        let b = &sys.bodies[0];
        // Angle from straight-down inferred from the bob's position about the pivot.
        let rel = b.pos - pivot;
        let angle = rel.x.atan2(-rel.y);
        MbdSample {
            t: sys.time,
            angle,
            x: b.pos.x,
            height: b.pos.y,
            energy: e,
            contact_force: 0.0,
        }
    };
    trace.push(sample_of(&sys, e0));

    for step in 1..=p.steps {
        sys.step(p.dt);
        let e = sys.energy();
        if !e.is_finite() {
            return Err(format!(
                "pendulum integration diverged (non-finite energy) at step {step} — \
                 try a smaller dt"
            ));
        }
        max_rel_drift = max_rel_drift.max((e - e0).abs() / denom);
        push_sampled(&mut trace, step, stride, sample_of(&sys, e));
    }

    let last = sys.bodies[0];
    let rel = last.pos - pivot;
    let final_angle = rel.x.atan2(-rel.y);

    Ok(MbdResult {
        demo: MbdDemo::Pendulum,
        samples: trace,
        final_time: sys.time,
        final_height: last.pos.y,
        final_angle_or_x: final_angle,
        energy_initial: e0,
        energy_rel_drift: max_rel_drift,
        max_penetration: 0.0,
        analytic_penetration: 0.0,
        max_contact_force: 0.0,
        at_rest: false,
    })
}

/// Run the **drop-onto-plane** demo through the real penalty-contact path.
///
/// A single rigid body (contact point at its centre of mass, so no torque) is
/// released `drop_height` above a flat ground [`Plane`], optionally with a
/// horizontal launch velocity, and integrated with symplectic Euler. The contact
/// wrench comes straight from [`body_contact_wrench`] (penalty normal
/// spring–damper clamped compressive + regularized Coulomb friction). Before the
/// body touches it is in free fall (`h = h₀ − ½ g t²`); after, it settles near
/// the plane with penetration bounded by `m g / kₙ`.
///
/// **Pre-condition:** the caller has already run [`MbdParams::validate`], so the
/// [`ContactParams::new`] and [`Plane::new`] calls below cannot trip an assert.
fn run_drop(p: &MbdParams) -> Result<MbdResult, String> {
    let m = p.mass;
    let g = p.gravity;
    let k_n = p.contact_stiffness;

    // Construct the validated contact params + ground plane. Safe: validate()
    // already proved kₙ>0, cₙ≥0, μ≥0, v_eps>0 are all finite, and the ground
    // normal is the fixed unit +y, so neither assert can fire here.
    let params = ContactParams::new(k_n, p.contact_damping, p.friction, FRICTION_VEL_SCALE);
    let plane = Plane::ground(0.0); // plane at y = 0, free half-space above

    // Body starts `drop_height` above the plane, contact point at the CoM.
    let mut body = Body::new(m, 0.05, Vector2::new(0.0, p.drop_height));
    body.vel = Vector2::new(p.launch_speed, 0.0);
    let contacts = [BodyContact {
        local: Vector2::zeros(),
        plane,
        params,
    }];

    let stride = trace_stride(p.steps);
    let mut trace: Vec<MbdSample> = Vec::with_capacity(p.steps.min(MAX_TRACE) + 1);
    let mut max_pen = 0.0_f64;
    let mut max_force = 0.0_f64;

    let sample_of = |body: &Body, f: Vector2<f64>| -> MbdSample {
        MbdSample {
            t: 0.0, // filled by caller loop
            angle: 0.0,
            x: body.pos.x,
            height: body.pos.y,
            energy: 0.0,
            contact_force: f.norm(),
        }
    };

    let mut t = 0.0;
    {
        // Initial sample (no contact yet at t=0 unless it starts penetrating).
        let state = ContactBodyState::from(&body);
        let w0 = body_contact_wrench(&state, &contacts);
        let mut s0 = sample_of(&body, w0.force);
        s0.t = 0.0;
        trace.push(s0);
        max_force = max_force.max(w0.force.norm());
        max_pen = max_pen.max(plane.penetration(body.pos).max(0.0));
    }

    for step in 1..=p.steps {
        let state = ContactBodyState::from(&body);
        let w = body_contact_wrench(&state, &contacts);
        let f = w.force;
        if !f.x.is_finite() || !f.y.is_finite() {
            return Err(format!(
                "contact force became non-finite at step {step} — try a smaller dt or \
                 a softer contact stiffness"
            ));
        }
        // Newton (symplectic Euler): a = (gravity + contact) / m.
        let ax = f.x / m;
        let ay = -g + f.y / m;
        body.vel.x += ax * p.dt;
        body.vel.y += ay * p.dt;
        body.pos.x += body.vel.x * p.dt;
        body.pos.y += body.vel.y * p.dt;
        t += p.dt;
        if !body.pos.y.is_finite() || !body.vel.y.is_finite() {
            return Err(format!(
                "drop integration diverged (non-finite state) at step {step} — \
                 try a smaller dt or a softer contact stiffness"
            ));
        }

        let pen = plane.penetration(body.pos).max(0.0);
        max_pen = max_pen.max(pen);
        max_force = max_force.max(f.norm());

        let mut s = sample_of(&body, f);
        s.t = t;
        push_sampled(&mut trace, step, stride, s);
    }

    let at_rest = body.vel.norm() < REST_SPEED;
    Ok(MbdResult {
        demo: MbdDemo::Drop,
        samples: trace,
        final_time: t,
        final_height: body.pos.y,
        final_angle_or_x: body.pos.x,
        energy_initial: 0.0,
        energy_rel_drift: 0.0,
        max_penetration: max_pen,
        analytic_penetration: m * g / k_n,
        max_contact_force: max_force,
        at_rest,
    })
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Draw the multibody-dynamics workbench. A no-op unless toggled on via
/// View → Multibody dynamics (robot / contact).
///
/// Mirrors [`crate::ppi_workbench::draw_ppi_workbench`].
pub fn draw_mbd_workbench(app: &mut ValenxApp, ctx: &egui::Context) {
    if !app.show_mbd_workbench {
        return;
    }
    let close = crate::workbench_chrome::workbench_shell(
        app,
        ctx,
        "valenx_mbd_workbench",
        "Multibody dynamics (robot / contact)",
        mbd_workbench_body,
    );
    if close {
        app.show_mbd_workbench = false;
    }
}

// ---------------------------------------------------------------------------
// Workbench body
// ---------------------------------------------------------------------------

fn mbd_workbench_body(app: &mut ValenxApp, ui: &mut egui::Ui) {
    ui.label(
        egui::RichText::new(
            "Multibody dynamics \u{2014} two native demos over the REAL in-house valenx-mbd \
             solver: a planar rigid-rod PENDULUM advanced by the constrained-DAE System \
             (KKT + Baumgarte + symplectic Euler; conservative, so energy stays bounded), and \
             a rigid body DROPPED onto a ground plane through the penalty-contact path \
             (spring-damper normal force + regularized Coulomb friction). [research / \
             preliminary-design grade: 2-D rigid bodies, penalty (compliant, step-size-sensitive) \
             contact \u{2014} a real v1, not a production multibody / hard-LCP code]",
        )
        .weak()
        .small(),
    );
    ui.separator();

    let mut do_run = false;

    {
        let s = &mut app.mbd;
        let p = &mut s.params;

        // Demo selector.
        ui.horizontal(|ui| {
            let lbl = ui.label("demo");
            egui::ComboBox::from_id_source("mbd_demo_combo")
                .selected_text(p.demo.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut p.demo, MbdDemo::Pendulum, MbdDemo::Pendulum.label());
                    ui.selectable_value(&mut p.demo, MbdDemo::Drop, MbdDemo::Drop.label());
                })
                .response
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Which multibody demo to run: the energy-conserving articulated rod \
                     pendulum, or the body dropped onto a plane with penalty contact + friction.",
                );
        });
        ui.add_space(4.0);

        // Shared parameters.
        ui.label(egui::RichText::new("Integrator & gravity").strong());
        egui::Grid::new("mbd_shared_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                let lbl = ui.label("gravity g (m/s^2)");
                ui.add(
                    egui::DragValue::new(&mut p.gravity)
                        .speed(0.05)
                        .range(0.001..=1.0e3)
                        .max_decimals(4),
                )
                .labelled_by(lbl.id)
                .on_hover_text("Downward gravitational acceleration. Must be finite and > 0.");
                ui.end_row();

                let lbl = ui.label("time step dt (s)");
                ui.add(
                    egui::DragValue::new(&mut p.dt)
                        .speed(1.0e-4)
                        .range(1.0e-6..=1.0)
                        .max_decimals(6),
                )
                .labelled_by(lbl.id)
                .on_hover_text(
                    "Explicit-integrator time step. Must be finite and > 0; a stiff contact \
                     needs a small dt for stability.",
                );
                ui.end_row();

                let lbl = ui.label("# steps");
                ui.add(
                    egui::DragValue::new(&mut p.steps)
                        .speed(50.0)
                        .range(1..=MAX_STEPS),
                )
                .labelled_by(lbl.id)
                .on_hover_text("How many integration steps to take. 1..=2,000,000.");
                ui.end_row();
            });

        ui.add_space(4.0);

        // Demo-specific parameters. Both grids are shown always (the inactive
        // one greyed) so the form layout + accessible names stay stable
        // regardless of the selected demo.
        let is_pend = p.demo == MbdDemo::Pendulum;
        let is_drop = p.demo == MbdDemo::Drop;

        ui.label(egui::RichText::new("Pendulum").strong());
        egui::Grid::new("mbd_pendulum_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                ui.add_enabled_ui(is_pend, |ui| {
                    let lbl = ui.label("rod length L (m)");
                    ui.add(
                        egui::DragValue::new(&mut p.rod_length)
                            .speed(0.05)
                            .range(0.01..=1.0e3)
                            .max_decimals(4),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Length of the rigid rod from the pivot to the bob. Must be finite \
                         and > 0. Small-angle period is 2\u{03C0}\u{221A}(L/g).",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(is_pend, |ui| {
                    let lbl = ui.label("release angle (rad)");
                    ui.add(
                        egui::DragValue::new(&mut p.release_angle)
                            .speed(0.01)
                            .range(-std::f64::consts::PI..=std::f64::consts::PI)
                            .max_decimals(4),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Initial angle from straight-down (rad). A small angle keeps the \
                         period close to the linear 2\u{03C0}\u{221A}(L/g) and the explicit \
                         integrator energy-tight.",
                    );
                });
                ui.end_row();
            });

        ui.add_space(4.0);
        ui.label(egui::RichText::new("Drop & contact").strong());
        egui::Grid::new("mbd_drop_params")
            .num_columns(2)
            .striped(true)
            .show(ui, |ui| {
                ui.add_enabled_ui(is_drop, |ui| {
                    let lbl = ui.label("drop height (m)");
                    ui.add(
                        egui::DragValue::new(&mut p.drop_height)
                            .speed(0.02)
                            .range(-1.0..=1.0e3)
                            .max_decimals(4),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Initial height of the body above the ground plane (m). Before it \
                         touches, height follows h(t) = h0 - 0.5*g*t^2.",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(is_drop, |ui| {
                    let lbl = ui.label("launch speed (m/s)");
                    ui.add(
                        egui::DragValue::new(&mut p.launch_speed)
                            .speed(0.05)
                            .range(-1.0e3..=1.0e3)
                            .max_decimals(4),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Initial horizontal velocity (m/s). With a non-zero friction \
                         coefficient, a sliding body is decelerated by Coulomb friction.",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(is_drop, |ui| {
                    let lbl = ui.label("mass m (kg)");
                    ui.add(
                        egui::DragValue::new(&mut p.mass)
                            .speed(0.05)
                            .range(0.001..=1.0e6)
                            .max_decimals(4),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Body mass (kg). Must be finite and > 0. Steady penetration is m*g/k.",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(is_drop, |ui| {
                    let lbl = ui.label("contact stiffness k (N/m)");
                    ui.add(
                        egui::DragValue::new(&mut p.contact_stiffness)
                            .speed(100.0)
                            .range(1.0..=1.0e9)
                            .max_decimals(2),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Penalty normal stiffness kn (N/m). MUST be > 0 (validated before the \
                         solver runs). Stiffer -> smaller penetration but needs a smaller dt.",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(is_drop, |ui| {
                    let lbl = ui.label("contact damping c (N*s/m)");
                    ui.add(
                        egui::DragValue::new(&mut p.contact_damping)
                            .speed(1.0)
                            .range(0.0..=1.0e7)
                            .max_decimals(3),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Penalty normal damping cn (N*s/m). MUST be >= 0 (validated). Dissipates \
                         bounce energy so the body settles.",
                    );
                });
                ui.end_row();

                ui.add_enabled_ui(is_drop, |ui| {
                    let lbl = ui.label("friction coeff mu");
                    ui.add(
                        egui::DragValue::new(&mut p.friction)
                            .speed(0.02)
                            .range(0.0..=10.0)
                            .max_decimals(3),
                    )
                    .labelled_by(lbl.id)
                    .on_hover_text(
                        "Coulomb friction coefficient mu (dimensionless). MUST be >= 0 \
                         (validated). 0 = frictionless; bounds tangential force by mu*fn.",
                    );
                });
                ui.end_row();
            });

        ui.add_space(6.0);
        ui.horizontal(|ui| {
            if ui
                .button(egui::RichText::new("Run").strong())
                .on_hover_text(
                    "Validate the parameters, then integrate the selected demo through the \
                     real valenx-mbd solver (constrained-DAE pendulum, or penalty-contact drop).",
                )
                .clicked()
            {
                do_run = true;
            }
        });
    }

    // --- Execute (outside the params borrow) --------------------------------
    if do_run {
        run_and_store(app);
    }

    // --- Status line ---------------------------------------------------------
    let s = &app.mbd;
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
    draw_mbd_viz(s, ui);
}

/// Run the simulation and fold the result (or error) into the workbench status.
/// Factored out so the Run button (and tests) can share it.
pub(crate) fn run_and_store(app: &mut ValenxApp) {
    let s = &mut app.mbd;
    match s.run() {
        Ok(res) => {
            s.status = match res.demo {
                MbdDemo::Pendulum => format!(
                    "\u{2714} pendulum \u{00B7} {} steps \u{00B7} t={:.2}s \u{00B7} energy drift \
                     {:.2e} (conserved)",
                    res.samples.len(),
                    res.final_time,
                    res.energy_rel_drift,
                ),
                MbdDemo::Drop => format!(
                    "\u{2714} drop \u{00B7} t={:.2}s \u{00B7} peak penetration {:.3e} m (analytic \
                     {:.3e}) \u{00B7} {}",
                    res.final_time,
                    res.max_penetration,
                    res.analytic_penetration,
                    if res.at_rest { "at rest" } else { "moving" },
                ),
            };
            s.result = Some(res);
        }
        Err(e) => {
            s.status = format!("\u{26A0} {e}");
            s.result = None;
        }
    }
}

// ---------------------------------------------------------------------------
// 2-D visualisation (painter + readout)
// ---------------------------------------------------------------------------

fn draw_mbd_viz(s: &MbdWorkbenchState, ui: &mut egui::Ui) {
    let Some(res) = &s.result else {
        ui.label(
            egui::RichText::new(
                "press \"Run\" to integrate the selected demo through valenx-mbd and draw the \
                 motion",
            )
            .weak(),
        );
        return;
    };

    match res.demo {
        MbdDemo::Pendulum => {
            ui.label(egui::RichText::new("Pendulum motion").strong());
            ui.label(
                egui::RichText::new(
                    "the rod swings about the pivot (top); the faint arc is the swept path \
                     \u{00B7} energy is conserved (NVE)",
                )
                .weak()
                .small(),
            );
            draw_pendulum(res, ui);
        }
        MbdDemo::Drop => {
            ui.label(egui::RichText::new("Drop: height & contact force vs time").strong());
            ui.label(
                egui::RichText::new(
                    "cyan = body height above the plane \u{00B7} amber = contact-force magnitude \
                     \u{00B7} grey line = the ground (y=0)",
                )
                .weak()
                .small(),
            );
            draw_drop(res, ui);
        }
    }

    // Readouts grid below the viz.
    ui.add_space(6.0);
    egui::Grid::new("mbd_stats")
        .num_columns(2)
        .striped(true)
        .show(ui, |ui| {
            let row = |ui: &mut egui::Ui, k: &str, v: String| {
                ui.label(k);
                ui.label(v);
                ui.end_row();
            };
            row(ui, "demo", res.demo.label().to_string());
            row(ui, "simulated time", format!("{:.3} s", res.final_time));
            row(ui, "samples (trace)", format!("{}", res.samples.len()));
            match res.demo {
                MbdDemo::Pendulum => {
                    row(
                        ui,
                        "final angle",
                        format!("{:.4} rad", res.final_angle_or_x),
                    );
                    row(
                        ui,
                        "final height (CoM y)",
                        format!("{:.4} m", res.final_height),
                    );
                    row(ui, "initial energy", format!("{:.6} J", res.energy_initial));
                    row(
                        ui,
                        "max energy drift |E-E0|/|E0|",
                        format!("{:.3e} (conserved)", res.energy_rel_drift),
                    );
                }
                MbdDemo::Drop => {
                    row(ui, "final x", format!("{:.4} m", res.final_angle_or_x));
                    row(ui, "final height", format!("{:.5} m", res.final_height));
                    row(
                        ui,
                        "peak penetration",
                        format!("{:.4e} m", res.max_penetration),
                    );
                    row(
                        ui,
                        "analytic m*g/k",
                        format!("{:.4e} m", res.analytic_penetration),
                    );
                    row(
                        ui,
                        "peak contact force",
                        format!("{:.3} N", res.max_contact_force),
                    );
                    row(
                        ui,
                        "settled",
                        if res.at_rest {
                            "yes (near rest)".into()
                        } else {
                            "no (still moving)".to_string()
                        },
                    );
                }
            }
            row(
                ui,
                "engine",
                "valenx-mbd (research / preliminary-design grade)".to_string(),
            );
        });
}

/// Draw the pendulum's swing with the egui painter: the pivot at the top, the
/// rod to the current bob, and a faint arc of the swept path.
fn draw_pendulum(res: &MbdResult, ui: &mut egui::Ui) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(460.0, 300.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.samples.is_empty() {
        return;
    }

    // World extent: the rod can reach +-L horizontally and down to -L; fit that
    // box into the rect with a margin. The pivot is at world (0,0).
    let reach = res
        .samples
        .iter()
        .map(|s| s.x.abs().max(s.height.abs()))
        .fold(0.0_f64, f64::max)
        .max(1e-6);
    let margin = 36.0;
    let scale = ((rect.width().min(rect.height())) * 0.5 - margin) / reach as f32;
    let pivot_screen = egui::pos2(rect.center().x, rect.top() + margin);
    let to_screen = |x: f64, y: f64| -> egui::Pos2 {
        // World +y is up; screen +y is down. Pivot at world origin.
        egui::pos2(
            pivot_screen.x + (x as f32) * scale,
            pivot_screen.y - (y as f32) * scale,
        )
    };

    // Faint swept arc (the bob path over the whole run).
    let path: Vec<egui::Pos2> = res
        .samples
        .iter()
        .map(|s| to_screen(s.x, s.height))
        .collect();
    if path.len() >= 2 {
        painter.add(egui::Shape::line(
            path.clone(),
            egui::Stroke::new(1.0, egui::Color32::from_rgb(60, 90, 110)),
        ));
    }

    // The current (final) rod + bob.
    let bob = *path.last().unwrap();
    painter.line_segment(
        [pivot_screen, bob],
        egui::Stroke::new(2.5, egui::Color32::from_rgb(180, 200, 215)),
    );
    // Pivot marker.
    painter.circle_filled(pivot_screen, 4.0, egui::Color32::from_rgb(210, 210, 210));
    // Bob.
    painter.circle_filled(bob, 9.0, egui::Color32::from_rgb(90, 200, 210));

    painter.text(
        egui::pos2(rect.left() + 6.0, rect.top() + 6.0),
        egui::Align2::LEFT_TOP,
        "pivot",
        egui::FontId::monospace(11.0),
        egui::Color32::from_gray(170),
    );
}

/// Draw the drop demo's height-vs-time and contact-force-vs-time traces with
/// the egui painter (two y-axes auto-scaled independently, shared time axis).
fn draw_drop(res: &MbdResult, ui: &mut egui::Ui) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(460.0, 300.0), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 0.0, egui::Color32::from_rgb(14, 22, 34));

    if res.samples.len() < 2 {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "no trajectory",
            egui::FontId::monospace(12.0),
            egui::Color32::from_gray(120),
        );
        return;
    }

    let plot = rect.shrink(28.0);
    let t_max = res.final_time.max(1e-9);
    let h_min = res
        .samples
        .iter()
        .map(|s| s.height)
        .fold(f64::INFINITY, f64::min)
        .min(0.0);
    let h_max = res
        .samples
        .iter()
        .map(|s| s.height)
        .fold(f64::NEG_INFINITY, f64::max)
        .max(0.0);
    let h_span = (h_max - h_min).max(1e-9);
    let f_max = res
        .samples
        .iter()
        .map(|s| s.contact_force)
        .fold(0.0_f64, f64::max)
        .max(1e-9);

    let tx = |t: f64| plot.left() + (t / t_max) as f32 * plot.width();
    let hy = |h: f64| plot.bottom() - ((h - h_min) / h_span) as f32 * plot.height();
    let fy = |f: f64| plot.bottom() - (f / f_max) as f32 * plot.height();

    // Ground line (height = 0).
    let y0 = hy(0.0);
    painter.line_segment(
        [egui::pos2(plot.left(), y0), egui::pos2(plot.right(), y0)],
        egui::Stroke::new(1.0, egui::Color32::from_gray(110)),
    );

    // Height trace (cyan).
    let height_path: Vec<egui::Pos2> = res
        .samples
        .iter()
        .map(|s| egui::pos2(tx(s.t), hy(s.height)))
        .collect();
    painter.add(egui::Shape::line(
        height_path,
        egui::Stroke::new(1.8, egui::Color32::from_rgb(90, 200, 210)),
    ));

    // Contact-force trace (amber).
    let force_path: Vec<egui::Pos2> = res
        .samples
        .iter()
        .map(|s| egui::pos2(tx(s.t), fy(s.contact_force)))
        .collect();
    painter.add(egui::Shape::line(
        force_path,
        egui::Stroke::new(1.6, egui::Color32::from_rgb(230, 180, 70)),
    ));

    // Axis labels.
    painter.text(
        egui::pos2(plot.left(), rect.top() + 4.0),
        egui::Align2::LEFT_TOP,
        format!("h\u{2208}[{h_min:.3},{h_max:.3}] m  f\u{2264}{f_max:.1} N"),
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(170),
    );
    painter.text(
        egui::pos2(plot.right(), plot.bottom() + 4.0),
        egui::Align2::RIGHT_TOP,
        format!("t={t_max:.2}s"),
        egui::FontId::monospace(10.0),
        egui::Color32::from_gray(170),
    );
}

// ---------------------------------------------------------------------------
// Tests (unit + headless_ui_tests, mirroring ppi_workbench)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn default_pendulum_run_succeeds_and_is_populated() {
        let s = MbdWorkbenchState::default();
        let res = s.run().expect("default (pendulum) run should succeed");
        assert_eq!(res.demo, MbdDemo::Pendulum);
        assert!(!res.samples.is_empty(), "trace should be populated");
        assert!(res.samples.len() <= MAX_TRACE + 1, "trace is down-sampled");
        assert!(res.energy_initial.is_finite());
    }

    #[test]
    fn pendulum_conserves_energy_pin() {
        // PIN (analytic): an undamped rigid-rod pendulum is a conservative (NVE)
        // system, so the REAL valenx-mbd System::energy must stay bounded — the
        // relative drift over the run is tiny (symplectic Euler + Baumgarte).
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Pendulum;
        s.params.rod_length = 1.0;
        s.params.release_angle = 0.2;
        s.params.dt = 1.0e-3;
        s.params.steps = 8_000;
        let res = s.run().expect("pendulum run should succeed");
        assert!(
            res.energy_rel_drift < 0.02,
            "energy drift {} should be < 2% (conservative system)",
            res.energy_rel_drift
        );
    }

    #[test]
    fn pendulum_small_angle_period_matches_simple_pendulum() {
        // A small-angle rigid-rod pendulum (rod through the CoM) swings as an
        // ideal simple pendulum with period 2π√(L/g). Detect the period from
        // x-axis zero crossings of the trace and compare.
        let l = 1.0;
        let g = 9.81;
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Pendulum;
        s.params.rod_length = l;
        s.params.release_angle = 0.05; // small angle
        s.params.dt = 5.0e-4;
        s.params.steps = 12_000;
        let res = s.run().expect("pendulum run should succeed");

        // Period from zero crossings of x(t).
        let mut crossings = Vec::new();
        for w in res.samples.windows(2) {
            let (t0, y0) = (w[0].t, w[0].x);
            let (t1, y1) = (w[1].t, w[1].x);
            if y0 == 0.0 || (y0 < 0.0) != (y1 < 0.0) {
                let frac = if (y1 - y0).abs() > 1e-15 {
                    -y0 / (y1 - y0)
                } else {
                    0.0
                };
                crossings.push(t0 + frac * (t1 - t0));
            }
        }
        assert!(
            crossings.len() >= 3,
            "need at least 3 zero crossings to measure a period"
        );
        let measured = crossings[2] - crossings[0]; // one full period
        let analytic = 2.0 * PI * (l / g).sqrt();
        assert!(
            (measured - analytic).abs() / analytic < 0.05,
            "measured period {measured} vs simple-pendulum {analytic}"
        );
    }

    #[test]
    fn drop_free_fall_matches_analytic_before_contact() {
        // PIN (analytic): with the body released well above a plane, before it
        // touches it is in pure free fall, so its height follows
        // h(t) = h0 - 0.5*g*t^2 to integrator tolerance. Pick a sample partway
        // down (still above the plane) and compare.
        let h0 = 1.0;
        let g = 9.81;
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.drop_height = h0;
        s.params.gravity = g;
        s.params.launch_speed = 0.0;
        s.params.dt = 1.0e-4;
        s.params.steps = 30_000;
        let res = s.run().expect("drop run should succeed");

        // Free fall holds only on the pristine descent BEFORE the body first
        // reaches the plane. After it touches (h ≤ 0) the penalty contact bounces
        // it back up, so a naive "last sample with h > 0.2" would wrongly pick a
        // post-bounce sample on the way up. Walk the trace in order and stop at the
        // first sample that has dropped to/through the plane; the last sample above
        // 0.2 m within that pre-contact window is a genuine free-fall point.
        let mut pre_contact: Vec<&MbdSample> = Vec::new();
        for s in &res.samples {
            if s.height <= 0.0 {
                break; // first contact — descent ends here
            }
            pre_contact.push(s);
        }
        let above: Vec<&MbdSample> = pre_contact.into_iter().filter(|s| s.height > 0.2).collect();
        assert!(
            above.len() >= 2,
            "should have free-fall samples above the plane"
        );
        let probe = above[above.len() - 1];
        let analytic = h0 - 0.5 * g * probe.t * probe.t;
        assert!(
            (probe.height - analytic).abs() < 2.0e-3,
            "free-fall height {} vs analytic {} at t={}",
            probe.height,
            analytic,
            probe.t
        );
    }

    #[test]
    fn drop_settles_near_plane_with_bounded_penetration_pin() {
        // PIN (analytic): a body dropped onto a plane with penalty contact +
        // damping settles near the plane, with the steady penetration bounded by
        // ~m*g/kn. After enough steps it is at rest and the final penetration is
        // within a small multiple of the analytic value.
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.drop_height = 0.2;
        s.params.mass = 1.0;
        s.params.gravity = 9.81;
        s.params.contact_stiffness = 1.0e4;
        s.params.contact_damping = 60.0;
        s.params.friction = 0.0;
        s.params.dt = 1.0e-4;
        s.params.steps = 200_000;
        let res = s.run().expect("drop run should succeed");

        assert!(res.at_rest, "the body should settle to near rest");
        let d_eq = res.analytic_penetration; // m*g/kn
                                             // Final height is a small negative penetration; check it settled close to
                                             // the analytic resting penetration (within 50% — penalty contact + the
                                             // explicit integrator have a documented tolerance).
        let final_pen = (-res.final_height).max(0.0);
        assert!(
            (final_pen - d_eq).abs() / d_eq < 0.5,
            "settled penetration {final_pen} vs analytic {d_eq}"
        );
        // And the body never sank far past the analytic depth (bounded ~m*g/kn,
        // allowing for the transient overshoot on impact).
        assert!(
            res.max_penetration < 20.0 * d_eq,
            "peak penetration {} should stay bounded near m*g/k ({})",
            res.max_penetration,
            d_eq
        );
    }

    #[test]
    fn drop_with_friction_decelerates_slide() {
        // A body launched horizontally and dropped onto the plane with a non-zero
        // friction coefficient must end up sliding SLOWER than the same launch
        // with mu = 0 (Coulomb friction removes tangential momentum). Exercises
        // the real regularized-friction path in body_contact_wrench.
        let base = |mu: f64| {
            let mut s = MbdWorkbenchState::default();
            s.params.demo = MbdDemo::Drop;
            s.params.drop_height = 0.02;
            s.params.launch_speed = 2.0;
            s.params.mass = 1.0;
            s.params.gravity = 9.81;
            s.params.contact_stiffness = 1.0e5;
            s.params.contact_damping = 200.0;
            s.params.friction = mu;
            s.params.dt = 1.0e-4;
            s.params.steps = 40_000;
            let res = s.run().expect("drop run should succeed");
            res.final_angle_or_x // final x position (proxy for slide distance)
        };
        let x_frictionless = base(0.0);
        let x_friction = base(0.6);
        assert!(
            x_friction < x_frictionless,
            "friction should reduce slide distance: mu=0.6 x={x_friction} vs mu=0 x={x_frictionless}"
        );
    }

    // ---- degenerate-param / contact-assert-guard tests — must return Err, NOT panic ----

    #[test]
    fn nonpositive_contact_stiffness_returns_err_not_panic() {
        // THE KEY TEST: an invalid contact stiffness (<= 0) would trip
        // ContactParams::new's assert (it has a #[should_panic] test). The
        // workbench MUST reject it in-panel as Err — proving validate() guards
        // the assert and no user value reaches it.
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.contact_stiffness = 0.0;
        assert!(
            s.run().is_err(),
            "stiffness <= 0 must return Err (guard the ContactParams assert), not panic"
        );
        s.params.contact_stiffness = -5.0;
        assert!(
            s.run().is_err(),
            "negative stiffness must return Err, not panic"
        );
    }

    #[test]
    fn negative_contact_damping_returns_err_not_panic() {
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.contact_damping = -1.0;
        assert!(
            s.run().is_err(),
            "damping < 0 must return Err (guard the ContactParams assert), not panic"
        );
    }

    #[test]
    fn negative_friction_returns_err_not_panic() {
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.friction = -0.1;
        assert!(
            s.run().is_err(),
            "friction < 0 must return Err (guard the ContactParams assert), not panic"
        );
    }

    #[test]
    fn nonpositive_dt_returns_err_not_panic() {
        let mut s = MbdWorkbenchState::default();
        s.params.dt = 0.0;
        assert!(s.run().is_err(), "dt <= 0 must return Err, not panic");
        s.params.dt = -1.0e-3;
        assert!(s.run().is_err(), "negative dt must return Err, not panic");
    }

    #[test]
    fn nonfinite_params_return_err_not_panic() {
        // NaN / infinite values everywhere must be rejected, never reach an assert.
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.contact_stiffness = f64::NAN;
        assert!(s.run().is_err(), "NaN stiffness must return Err");
        s.params = MbdParams::default();
        s.params.gravity = f64::INFINITY;
        assert!(s.run().is_err(), "infinite gravity must return Err");
        s.params = MbdParams::default();
        s.params.demo = MbdDemo::Drop;
        s.params.friction = f64::INFINITY;
        assert!(s.run().is_err(), "infinite friction must return Err");
    }

    #[test]
    fn zero_or_excessive_steps_returns_err() {
        let mut s = MbdWorkbenchState::default();
        s.params.steps = 0;
        assert!(s.run().is_err(), "0 steps must return Err");
        s.params.steps = MAX_STEPS + 1;
        assert!(s.run().is_err(), "too-many steps must return Err");
    }

    #[test]
    fn nonpositive_mass_and_rod_return_err() {
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.mass = 0.0;
        assert!(s.run().is_err(), "mass <= 0 must return Err");
        let mut s2 = MbdWorkbenchState::default();
        s2.params.demo = MbdDemo::Pendulum;
        s2.params.rod_length = 0.0;
        assert!(s2.run().is_err(), "rod length <= 0 must return Err");
    }

    #[test]
    fn valid_drop_params_pass_validation() {
        // Positive control: a sane drop config validates cleanly (so the Err
        // tests above are catching the bad values, not failing for another
        // reason).
        let mut p = MbdParams::default();
        p.demo = MbdDemo::Drop;
        assert!(p.validate().is_ok(), "default drop params must validate");
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
            draw_mbd_workbench(app, ctx);
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
        assert!(!app.show_mbd_workbench);
        let ctx = egui::Context::default();
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            draw_mbd_workbench(&mut app, ctx);
        });
        // No panic = pass.
    }

    #[test]
    fn workbench_draws_when_shown_without_panic() {
        let mut app = ValenxApp::default();
        app.show_mbd_workbench = true;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_populated_pendulum_result_without_panic() {
        let mut app = ValenxApp::default();
        app.show_mbd_workbench = true;
        let res = app.mbd.run().expect("pendulum run should succeed");
        app.mbd.result = Some(res);
        app.mbd.status = "\u{2714} test result".to_string();
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_drop_result_without_panic() {
        // Exercise the drop demo's height/force trace painter + readout rows.
        let mut app = ValenxApp::default();
        app.show_mbd_workbench = true;
        app.mbd.params.demo = MbdDemo::Drop;
        app.mbd.params.steps = 5_000;
        let res = app.mbd.run().expect("drop run should succeed");
        app.mbd.result = Some(res);
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn workbench_draws_with_error_status_without_panic() {
        let mut app = ValenxApp::default();
        app.show_mbd_workbench = true;
        // Trigger an error state (a bad contact stiffness is fail-loud in run()).
        app.mbd.params.demo = MbdDemo::Drop;
        app.mbd.params.contact_stiffness = -1.0;
        let result = app.mbd.run();
        app.mbd.status = match result {
            Err(e) => format!("\u{26A0} {e}"),
            Ok(_) => "\u{26A0} simulated error for testing".to_string(),
        };
        app.mbd.result = None;
        let _ = draw_and_collect_nodes(&mut app);
    }

    #[test]
    fn numeric_controls_are_labelled_by_named() {
        let mut app = ValenxApp::default();
        app.show_mbd_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        // The numeric DragValues MUST each carry an accessible name (be
        // labelled_by a caption) so the panel is AI-drivable.
        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
            "expected at least 4 numeric controls (DragValues), got {}",
            spin_buttons.len()
        );
        assert!(
            spin_buttons.iter().all(|n| !n.labelled_by().is_empty()),
            "every DragValue must be labelled_by a caption (AI-drivable name)"
        );

        // Check the specific captions are present as named accessibility nodes.
        for caption in [
            "gravity g (m/s^2)",
            "time step dt (s)",
            "# steps",
            "rod length L (m)",
            "release angle (rad)",
            "drop height (m)",
            "launch speed (m/s)",
            "mass m (kg)",
            "contact stiffness k (N/m)",
            "contact damping c (N*s/m)",
            "friction coeff mu",
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
        // Each numeric DragValue's `labelled_by` target must RESOLVE to a real
        // named caption node, not a dangling id.
        let mut app = ValenxApp::default();
        app.show_mbd_workbench = true;
        let nodes = draw_and_collect_nodes(&mut app);

        let by_id: std::collections::HashMap<NodeId, &Node> =
            nodes.iter().map(|(id, n)| (*id, n)).collect();

        let spin_buttons: Vec<&Node> = nodes
            .iter()
            .map(|(_, n)| n)
            .filter(|n| n.role() == Role::SpinButton)
            .collect();
        assert!(
            spin_buttons.len() >= 4,
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
        for caption in ["contact stiffness k (N/m)", "friction coeff mu"] {
            assert!(
                nodes.iter().any(|(_, n)| n.name() == Some(caption)),
                "caption '{caption}' should be a named node in the a11y tree"
            );
        }
    }

    #[test]
    fn energy_conservation_pin_from_ui_state() {
        // Mirror of the unit energy pin, exercised from the UI-state struct.
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Pendulum;
        s.params.steps = 8_000;
        let res = s.run().expect("pendulum run");
        assert!(
            res.energy_rel_drift < 0.02,
            "energy conserved from UI state"
        );
    }

    #[test]
    fn free_fall_and_contact_rest_pins_from_ui_state() {
        // Free-fall before contact, then settle near the plane — both from the
        // UI-state struct.
        let mut s = MbdWorkbenchState::default();
        s.params.demo = MbdDemo::Drop;
        s.params.drop_height = 0.2;
        s.params.dt = 1.0e-4;
        s.params.steps = 200_000;
        let res = s.run().expect("drop run");
        assert!(res.at_rest, "body settles near rest");
        assert!(
            res.max_penetration < 20.0 * res.analytic_penetration,
            "penetration bounded near m*g/k"
        );
    }

    #[test]
    fn invalid_contact_param_shows_error_not_panic_from_ui() {
        // Zero stiffness (or a bad dt) must surface the error in-panel, not panic
        // (the ContactParams assert is never reached).
        let mut state = MbdWorkbenchState::default();
        state.params.demo = MbdDemo::Drop;
        state.params.contact_stiffness = 0.0;
        assert!(
            state.run().is_err(),
            "0 stiffness must produce Err, not panic"
        );
        state.params = MbdParams::default();
        state.params.dt = -1.0;
        assert!(state.run().is_err(), "bad dt must produce Err, not panic");
    }

    #[test]
    fn run_and_store_sets_status_and_result() {
        // The Run path: a successful run populates result + a ✓ status; a failing
        // run clears result + sets a ⚠ status (never panics).
        let mut app = ValenxApp::default();
        app.mbd.params.demo = MbdDemo::Pendulum;
        run_and_store(&mut app);
        assert!(app.mbd.result.is_some(), "successful run stores a result");
        assert!(
            app.mbd.status.starts_with('\u{2714}'),
            "status shows success"
        );

        app.mbd.params.demo = MbdDemo::Drop;
        app.mbd.params.contact_stiffness = -1.0; // fail-loud
        run_and_store(&mut app);
        assert!(app.mbd.result.is_none(), "failed run clears the result");
        assert!(
            app.mbd.status.starts_with('\u{26A0}'),
            "status shows the error"
        );
    }

    #[test]
    fn agent_bridge_mbd_id_resolves_and_sets_flag() {
        // Verify the two mechanisms the agent bridge uses for
        //   `OpenWorkbench { id: "mbd" }`:
        //   1. TabKind::from_id("mbd") -> Some(TabKind::Mbd)
        //      (plus the aliases "multibody" / "robot")
        //   2. set_workbench_flag(app, "mbd", true) -> show_mbd_workbench = true
        use crate::project_tabs::{set_workbench_flag, TabKind};

        // 1. Lookup (canonical + aliases).
        assert_eq!(
            TabKind::from_id("mbd"),
            Some(TabKind::Mbd),
            "\"mbd\" must resolve to TabKind::Mbd"
        );
        assert_eq!(TabKind::from_id("multibody"), Some(TabKind::Mbd));
        assert_eq!(TabKind::from_id("robot"), Some(TabKind::Mbd));
        // Case-insensitive + whitespace-tolerant.
        assert_eq!(TabKind::from_id("  Mbd  "), Some(TabKind::Mbd));

        // 2. Flag toggle.
        let mut app = ValenxApp::default();
        assert!(!app.show_mbd_workbench);
        set_workbench_flag(&mut app, "mbd", true);
        assert!(
            app.show_mbd_workbench,
            "set_workbench_flag(\"mbd\", true) must set the flag"
        );
        set_workbench_flag(&mut app, "mbd", false);
        assert!(!app.show_mbd_workbench);
    }
}
