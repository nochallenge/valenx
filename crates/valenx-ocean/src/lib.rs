//! # valenx-ocean — ocean wave field and rigid-body buoyancy
//!
//! An **in-house, dependency-light** library for an ocean **wave field** and the
//! **quasi-static buoyancy** of a floating rigid body. Two pieces:
//!
//! * [`wave`] — a **sum-of-`N` directional Gerstner (trochoidal) wave field**.
//!   Each component has a wavelength, amplitude, direction, steepness, and phase
//!   and obeys the **deep-water dispersion relation** `ω = sqrt(g k)` with
//!   `k = 2π / L`. The field reports surface **displacement** (the 3-D offset of
//!   a rest point), **water height**, and the analytic **surface normal** at a
//!   horizontal `(x, z)` point and time `t`. (An FFT/Tessendorf spectral ocean
//!   — a Phillips/JONSWAP-seeded inverse-FFT height field — is a documented
//!   follow-up.)
//! * [`buoyancy`] — **Archimedes buoyancy** on a floating body. A body is either
//!   a set of [sample volume points](buoyancy::SampleBody) or a simple analytic
//!   [prismatic hull](buoyancy::HullBody) (waterplane area + volume). Against the
//!   wave surface the model computes the **submerged volume**, the **buoyant
//!   force** `ρ_water · g · V_submerged` at the **centre of buoyancy**, a
//!   **linear + quadratic drag**, the **net force/torque**, and integrates a
//!   **heave/pitch/roll** response.
//!
//! ## Clean-room provenance
//!
//! This is a **clean-room reimplementation from the published math** — Gerstner
//! / trochoidal deep-water waves (a classical result; see Tessendorf 2001,
//! *"Simulating Ocean Water"*, for the Gerstner formulation and surface
//! Jacobian) and Archimedes' principle. It is *algorithmically inspired* by the
//! existence of the **deprecated, Unreal-tied UE4 `OceanProject`** but is **not
//! a port** of it: that project is C++/Blueprint bound to the Unreal Engine and
//! not portable, so nothing was copied — the wave sum, dispersion relation,
//! surface normal, and buoyancy integrals are coded directly from the equations.
//!
//! ## ⚠ Honesty / scope — Gerstner + quasi-static buoyancy, NOT seakeeping CFD
//!
//! **This is a graphics / first-cut-engineering ocean: a Gerstner wave field with
//! quasi-static (hydrostatic) buoyancy. It is NOT a seakeeping CFD/RANS solver.**
//! In particular:
//!
//! * The waves are a **prescribed kinematic sum of trochoids**, not a solution of
//!   the free-surface Navier–Stokes (or even potential-flow) equations; there is
//!   no wave breaking, no wind input/dissipation balance, no nonlinear
//!   wave–wave energy transfer.
//! * Buoyancy is **Froude–Krylov / hydrostatic only**: the force comes from the
//!   undisturbed wave height under each sample point. There is **no diffraction**
//!   (the body does not scatter the incident wave) and **no radiation** (the
//!   body's motion does not generate outgoing waves).
//! * **Added mass and wave-radiation damping are not modelled** beyond a lumped
//!   linear+quadratic [`Drag`] coefficient — there is no
//!   frequency-dependent `A(ω)`/`B(ω)` from a boundary-element or strip-theory
//!   solve, and the rotational response uses a single scalar inertia.
//!
//! What *is* pinned in the tests are the **exact analytic checks** the model must
//! obey — not a comparison to a seakeeping reference.
//!
//! ## Pinned benchmarks (fail-loud, never fudged)
//!
//! 1. **Single Gerstner wave** *(`wave::tests`)* — crest-to-trough height is
//!    exactly `2·amplitude`; the deep-water phase speed is `c = sqrt(g / k)`; the
//!    surface is periodic in space (period `L`) and in time (period
//!    `T = 2π / ω`) to `< 1e-9`.
//! 2. **Archimedes** *(`buoyancy::tests`)* — a fully-submerged body of volume
//!    `V` feels `ρ g V` exactly; a body floating at equilibrium displaces its own
//!    weight (`V_sub · ρ · g = m · g`).
//! 3. **Heave natural frequency** *(`buoyancy::tests`)* — a prismatic body
//!    released from rest near equilibrium oscillates about its waterline at the
//!    analytic `ω_n = sqrt(ρ g A_wp / m)` (restoring `∝ ρ g A_wp`), within ~1 %.
//! 4. **Zero net force at equilibrium** *(`buoyancy::tests`)* — at the floating
//!    draft buoyancy cancels weight and the body does not drift.
//!
//! ## Fail-loud configuration & guarded divides
//!
//! Every constructor validates its parameters and returns
//! [`Result<_, OceanError>`] (never a panic or a silent `NaN`): a non-positive
//! wavelength, amplitude, steepness outside `[0, 1]`, a non-positive water
//! density, gravity, body mass, waterplane area, hull volume, inertia, or time
//! step; a negative drag coefficient; a zero-length wave direction; and any
//! non-finite (`NaN`/`±∞`) input. Every divide is guarded — the wavenumber
//! `k > 0`, the angular frequency `ω > 0`, the mass `m > 0`, the inertia `I > 0`,
//! and the density `ρ > 0` by construction.
//!
//! ## No `unsafe`, no `rand`
//!
//! Pure Rust with `#![forbid(unsafe_code)]`. Wave phases for the
//! [`deterministic_sea`](wave::OceanWaveField::deterministic_sea) preset come
//! from an **in-crate deterministic integer hash** (a SplitMix64 finaliser), so
//! results are reproducible on every run and platform without a `rand`
//! dependency.
//!
//! ## Connects to
//!
//! Marine / Navy hydro: a fast, deterministic wave-and-float model for vehicle
//! and platform behaviour at the sea surface — buoy/USV bobbing, a first-cut hull
//! response, sea-state visualisation — feeding the broader Valenx vehicle and
//! environment models.

#![forbid(unsafe_code)]

pub mod buoyancy;
pub mod error;
pub mod wave;

pub use buoyancy::{archimedes_force, BodyState, BuoyancySim, Drag, HullBody, Loads, SampleBody};
pub use error::OceanError;
pub use wave::{GerstnerWave, OceanWaveField, SEAWATER_DENSITY, STANDARD_GRAVITY};
