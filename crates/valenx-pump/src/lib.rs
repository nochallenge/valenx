//! # valenx-pump
//!
//! Closed-form **centrifugal-pump hydraulics** — the textbook handful of
//! relations that connect a pump's speed, the flow it delivers, the head
//! it produces, the power it draws, and whether it will cavitate.
//!
//! ## What this is
//!
//! A small, dependency-light library of the canonical centrifugal-pump
//! equations you find in any fluid-mechanics or turbomachinery text
//! (White, *Fluid Mechanics*; Munson; the Hydraulic Institute primers).
//! Everything is an explicit algebraic formula evaluated in SI units —
//! there is no CFD, no impeller geometry, no empirical performance map.
//!
//! Five topics, one module each:
//!
//! - **Affinity laws** ([`affinity`]) — how a fixed-geometry pump scales
//!   with shaft speed: flow `Q ∝ N`, head `H ∝ N²`, shaft power
//!   `P ∝ N³`. Run a measured duty point ([`affinity::DutyPoint`]) to a
//!   new speed with [`affinity::scale_to_speed`], and read off the
//!   speed-independent [`affinity::dimensionless_specific_speed`]
//!   `Ω_s = ω·√Q/(g·H)^¾` that classifies the impeller type.
//! - **System curve** ([`system`]) — the resistance the piping presents,
//!   `H = H_static + K·Q²`: a fixed static lift plus a velocity-head loss
//!   that grows with the square of flow ([`system::SystemCurve`]).
//! - **Hydraulic power** ([`power`]) — the useful power added to the
//!   fluid, `P = ρ·g·Q·H`, and the shaft power that follows once an
//!   efficiency is applied ([`power::hydraulic_power_w`],
//!   [`power::shaft_power_w`]).
//! - **NPSH** ([`npsh`]) — the available net positive suction head,
//!   `NPSHa = (P_atm − P_vap)/(ρ·g) + H_s − H_loss`, and the
//!   cavitation-margin check against a pump's required NPSH
//!   ([`npsh::SuctionConditions`], [`npsh::available_npsh_m`]).
//! - **Operating point** ([`operating`]) — where a pump's head curve and
//!   the system curve cross, found by closed-form intersection of a
//!   quadratic pump curve `H = H₀ − a·Q²` with the system curve
//!   ([`operating::operating_point`]).
//!
//! ```
//! use valenx_pump::{
//!     affinity::{scale_to_speed, DutyPoint},
//!     power::hydraulic_power_w,
//!     system::SystemCurve,
//! };
//!
//! // A pump tested at 1450 rpm delivers 0.05 m³/s at 30 m of head,
//! // drawing 18 kW. What does it do at 2900 rpm (double speed)?
//! let base = DutyPoint::new(1450.0, 0.05, 30.0, 18_000.0).unwrap();
//! let fast = scale_to_speed(&base, 2900.0).unwrap();
//! assert!((fast.flow_m3s - 0.10).abs() < 1e-12); //   Q doubles
//! assert!((fast.head_m - 120.0).abs() < 1e-9); //     H ×4
//! assert!((fast.power_w - 144_000.0).abs() < 1e-6); // P ×8
//!
//! // Useful power added to water at the fast point:
//! let p = hydraulic_power_w(1000.0, fast.flow_m3s, fast.head_m).unwrap();
//! assert!((p - 1000.0 * 9.80665 * 0.10 * 120.0).abs() < 1e-6);
//!
//! // A system that lifts 10 m and loses head as 2000·Q²:
//! let sys = SystemCurve::new(10.0, 2000.0).unwrap();
//! assert!((sys.head_m(0.10) - (10.0 + 2000.0 * 0.10 * 0.10)).abs() < 1e-12);
//! ```
//!
//! ## Model
//!
//! Conventions and the exact closed forms used throughout:
//!
//! - **Units are SI.** Flow `Q` in m³/s, head `H` in metres of the
//!   pumped fluid, power `P` in watts, pressures in pascals, density `ρ`
//!   in kg/m³, speed `N` in any consistent unit (rpm or rad/s — only the
//!   *ratio* enters the affinity laws). The gravitational acceleration is
//!   the standard [`G`] = 9.80665 m/s².
//! - **Affinity (fixed geometry, same fluid):**
//!   `Q₂/Q₁ = N₂/N₁`, `H₂/H₁ = (N₂/N₁)²`, `P₂/P₁ = (N₂/N₁)³`.
//! - **Specific speed:** `Ω_s = ω·√Q / (g·H)^(3/4)` with `ω` in rad/s —
//!   the dimensionless impeller-shape parameter, invariant under the
//!   affinity laws above (the `N·N^½/N^(3/2)` factors cancel).
//! - **System curve:** `H_sys(Q) = H_static + K·Q²` with resistance
//!   coefficient `K ≥ 0` (units m·s²/m⁶) and static head `H_static`
//!   (which may be negative for a flooded-suction / downhill system).
//! - **Hydraulic (water / fluid) power:** `P_hyd = ρ·g·Q·H`. Dividing by
//!   a pump efficiency `η ∈ (0, 1]` gives the shaft power
//!   `P_shaft = P_hyd / η`.
//! - **Available NPSH:** `NPSHa = (P_atm − P_vap)/(ρ·g) + H_s − H_loss`,
//!   where `H_s` is the static suction head (positive for flooded
//!   suction, **negative** for a suction lift) and `H_loss ≥ 0` is the
//!   friction loss in the suction line. The cavitation margin is
//!   `NPSHa − NPSHr`; a pump is safe when it is positive.
//! - **Operating point:** with the pump curve written as the downward
//!   parabola `H_pump(Q) = H₀ − a·Q²` (`a > 0`) and the system curve
//!   `H_sys(Q) = H_static + K·Q²`, equating heads gives the unique
//!   non-negative root
//!   `Q* = sqrt((H₀ − H_static)/(a + K))`, with
//!   `H* = H_static + K·Q*²`.
//!
//! ## Honest scope
//!
//! These are the genuine textbook relations — the affinity exponents are
//! exact, hydraulic power is exact, the system curve and the operating
//! point are solved in closed form, and every result is checked against
//! hand-computed ground truth in the test suite. The library is
//! deliberately a **closed-form teaching / preliminary-sizing tool**, and
//! it stops well short of a production pump-selection package:
//!
//! - **No real performance map.** A manufacturer's H–Q curve is empirical
//!   and droops in a way no single parabola captures; here the pump curve
//!   is the idealised `H₀ − a·Q²` parabola. Efficiency is a single user
//!   number, not an efficiency-island map, so the best-efficiency point
//!   and part-load penalties are not modelled.
//! - **Affinity laws assume ideal scaling** — fixed geometry, the same
//!   fluid, fully turbulent flow, and constant efficiency between speeds.
//!   Real Reynolds-number and slip effects make the cubed power law only
//!   approximate, especially over large speed ratios or with impeller
//!   trimming (which scales differently from speed).
//! - **NPSH and the system curve are lumped.** `H_loss` and `K` are
//!   single coefficients you supply; the crate does not compute pipe
//!   friction from Darcy–Weisbach, fittings, two-phase effects, or
//!   transient water-hammer. NPSHr is a pump datum you pass in, not a
//!   prediction.
//! - **Incompressible, single-phase, steady.** No slurry, gas
//!   entrainment, viscosity correction, multistage staging, or
//!   parallel/series pump combination logic.
//!
//! Treat the numbers as first-pass engineering estimates for learning and
//! preliminary sizing — **not** as a substitute for a validated selection
//! program, a vendor performance curve, or a hydraulic-transient study.
//!
//! Research/educational grade: textbook closed-form/numerical models; NOT
//! a clinical/medical/production engineering tool.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod affinity;
pub mod error;
pub mod npsh;
pub mod operating;
pub mod power;
pub mod system;

pub use error::PumpError;

/// Standard gravitational acceleration, in metres per second squared.
///
/// The conventional value adopted by the General Conference on Weights
/// and Measures (CGPM, 1901). Used wherever a head (in metres of fluid)
/// is converted to or from a pressure or a power.
pub const G: f64 = 9.806_65;
