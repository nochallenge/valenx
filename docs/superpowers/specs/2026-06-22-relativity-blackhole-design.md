# valenx-relativity — in-house general-relativity & black-hole engine

**Date:** 2026-06-22
**Lane:** Chat 2 (science/sim), branch `feat/science-lane`
**Status:** approved by user; build autonomously (full ladder, code-review until clean, ship via PRs).

## Identity

A new pure-Rust leaf crate `valenx-relativity`. It is a general-relativity
tensor-calculus engine; black holes are its built-in spacetimes and headline
showcase.

> In-house general-relativity engine: automatic-differentiation tensor calculus
> (Christoffel / Riemann / Ricci / Einstein, Kretschmann) over arbitrary metrics,
> the Schwarzschild / Kerr / Reissner–Nordström / Kerr–Newman black-hole family,
> geodesic integration, and black-hole shadow ray-tracing — every quantity
> cross-checked against closed-form ground truth.

In-house, not a wrapper: we do **not** shell out to SageMath / Mathematica /
Einstein Toolkit. The tensor calculus is built from scratch in Rust.

## Conventions

- **Geometrized units G = c = 1** internally. Mass is a length; the natural scale
  is `M = 1`, so radii come out in units of `M` (horizon at `r = 2`). A thin
  `units.rs` layer converts to/from SI (solar masses, kg, km, K, years).
- **Signature** (−,+,+,+). Coordinates: Boyer–Lindquist `(t, r, θ, φ)` for the
  Kerr–Newman family; spherical `(t, r, θ, φ)` for Schwarzschild/RN; Cartesian
  for Minkowski.
- **Honesty:** every fallible path returns `Result`; we never emit a silent
  NaN/Inf (matches the repo's astro-loop "seal silent NaN" hardening). Quantities
  with no closed form available are documented as such, not faked.

## Differentiation strategy (the core decision)

**Forward-mode automatic differentiation** (approach B), with hand-coded analytic
curvature used only as a *test oracle* (approach C).

- `Dual` carries value + first derivative; `HyperDual` carries value + two first
  derivatives + the mixed/second derivative. One HyperDual pass over the generic
  `metric<T>` yields `g`, `∂ₐg`, `∂ₐ∂_b g` exactly (machine precision).
- Christoffel symbols and the Riemann tensor are then **assembled from the
  standard formulas** using those partials — we do *not* differentiate through
  the matrix inverse (robust, textbook method).
- Finite differences (approach A) remain documented as a fallback.

## Module layout (== bottom-up PR ladder)

```
crates/valenx-relativity/
  src/
    lib.rs            crate docs, re-exports, RelativityError
    autodiff/         Dual (1st deriv) + HyperDual (2nd deriv), no external dep
    tensor.rs         4-vectors, symmetric 4×4 metric, raise/lower, inverse (nalgebra)
    metric.rs         Spacetime trait: metric<T>(x:[T;4]) -> Mat4<T>; coord tag
    spacetimes/       kerr_newman (master) + schwarzschild/kerr/reissner_nordstrom
                      constructors + minkowski (flat baseline)
    curvature.rs      Christoffel, Riemann, Ricci, Ricci scalar, Einstein, Kretschmann
    observables.rs    horizons, ergosphere, photon sphere, ISCO, redshift, shadow radius
    thermo.rs         surface gravity, Hawking temp, Bekenstein–Hawking entropy, evaporation
    geodesics.rs      timelike/null geodesic integrator (RKF45); deflection, precession, capture
    shadow.rs         backward ray-trace → shadow + photon-ring image data, redshift map
    units.rs          SI ↔ geometrized
  tests/ground_truth.rs   the validation table
  examples/shadow.rs      render a shadow to an image for demos
```

Deps: `nalgebra`, `thiserror`, `serde` (+ `serde_json` dev) — same as sibling
pure crates. Pure-Rust, no GUI/`rfd` → joins `scripts/qa.sh` `PURE_CRATES`,
validated headless via `cargo test -p valenx-relativity`.

`KerrNewman { mass, spin a, charge Q }` is the master metric; `schwarzschild`,
`kerr`, `reissner_nordstrom` are named constructors (Q=0, a=0, or both). One code
path, four validated special cases.

## Public API (sketch)

```rust
pub trait Spacetime {
    fn metric<T: Scalar>(&self, x: [T; 4]) -> Mat4<T>;  // symmetric g_μν
    fn coords(&self) -> CoordSystem;
    fn mass(&self) -> f64;
}
pub struct KerrNewman { pub mass: f64, pub spin: f64, pub charge: f64 }
pub fn schwarzschild(mass: f64) -> KerrNewman;
pub fn kerr(mass: f64, spin: f64) -> KerrNewman;
pub fn reissner_nordstrom(mass: f64, charge: f64) -> KerrNewman;
pub struct Minkowski;

pub fn curvature_at(s: &impl Spacetime, x: [f64;4]) -> Result<Curvature>;
// Curvature { christoffel, riemann, ricci, ricci_scalar, einstein, kretschmann }

pub fn horizons(s) -> Result<Horizons>;        // r±, ergosphere
pub fn photon_sphere(s) -> Result<f64>;
pub fn isco(s, sense: OrbitSense) -> Result<f64>;
pub fn gravitational_redshift(s, r_emit, r_obs) -> Result<f64>;
pub fn shadow_radius(s) -> Result<f64>;
pub fn thermo(s) -> Result<Thermo>;            // κ, T_H, S, evaporation
pub fn integrate_geodesic(s, init, kind, opts) -> Result<Trajectory>;
pub fn light_deflection(s, impact_b) -> Result<f64>;
pub fn perihelion_precession(s, a, e) -> Result<f64>;
pub fn render_shadow(s, camera, w, h) -> Result<ShadowImage>;
```

Result structs derive `Serialize` for the app handoff.

## Error model

`thiserror` enum `RelativityError`:
- `InvalidParameter` — mass ≤ 0; super-extremal `a² + Q² > M²` rejected by default
  (naked singularity), with explicit `allow_naked` opt-in.
- `CoordinateSingularity` — on horizon/axis where Boyer–Lindquist coords blow up.
- `GeodesicNonConvergence` — integrator failed; returns partial trajectory + reason.
- `OutsideDomain` — point not in the metric's valid region.

## Validation table (closed-form ground truth)

| Quantity | Ground truth (G=c=1) | Test |
|---|---|---|
| Minkowski curvature | exactly 0 | all components < 1e-12 |
| Schwarzschild Ricci | 0 (vacuum) | each `Rμν` < 1e-9 |
| Kerr Ricci | 0 (vacuum) | each `Rμν` < 1e-8 |
| Schwarzschild Kretschmann | `48 M²/r⁶` | rel. err < 1e-9 at several r |
| AD Christoffel/Riemann | = hand-coded analytic | < 1e-10 (Schwarzschild, Kerr) |
| Schwarzschild horizon | `2M` | exact |
| Photon sphere (Schw.) | `3M` | exact |
| ISCO (Schw.) | `6M` | exact |
| Shadow radius (Schw.) | `3√3 M ≈ 5.196 M` | rel. err < 1e-6 |
| Kerr horizons | `M ± √(M²−a²)` | exact |
| Kerr ISCO (Bardeen) | Bardeen 1972 formula | < 1e-6; extremal prograde → M |
| RN horizons | `M ± √(M²−Q²)` | exact; Q=0 → Schwarzschild |
| Kerr–Newman horizons | `M ± √(M²−a²−Q²)` | a=0 → RN, Q=0 → Kerr |
| Surface gravity (Schw.) | `κ = 1/(4M)` | exact |
| Hawking temp (Schw.) | `T = 1/(8πM)` → 6.17e-8 K·(M_⊙/M) | rel. err < 1e-9 |
| BH entropy (Schw.) | `S = A/4 = 4πM²` | exact |
| Light deflection (weak) | `4M/b` → 1.75″ grazing Sun | rel. err < 1e-3 |
| Mercury precession | `6πM/(a(1−e²))` → 43″/century | rel. err < 1e-2 |
| Photon capture | `b < 3√3 M` captured | boundary < 1e-3 |

## PR ladder (each green via `cargo test -p valenx-relativity`, clippy clean, code-reviewed)

1. **Math core** — autodiff, tensors, `Spacetime` trait, all 5 metrics, `curvature.rs`;
   ground-truth tests (flat, vacuum, Kretschmann, AD-vs-analytic). Register in
   workspace `Cargo.toml` + `qa.sh` `PURE_CRATES`.
2. **Observables + thermodynamics + units** — horizons/ergosphere/photon sphere/
   ISCO/redshift/shadow radius, `thermo.rs`, SI layer + tests.
3. **Geodesics** — RKF45 integrator, orbits, deflection (1.75″), precession
   (43″/century), capture.
4. **Shadow / imaging** — backward ray-tracer → shadow + photon ring (Schwarzschild
   circular, Kerr asymmetric), redshift map; `examples/shadow.rs` PNG demo.

After PR4 (or as layers land), hand off to Chat 1 via `valenx-COORDINATION.md` for
the valenx-app workbench wiring (their lane).

## Out of scope (YAGNI for v1)

- Numerical relativity / dynamical spacetime evolution (BSSN, binary inspiral).
- Gravitational-wave waveforms / quasi-normal modes.
- A symbolic computer-algebra system (we use AD, not symbolic).
- Greybody factors / full Hawking spectrum (evaporation time uses the standard
  Page-style `∝ M³` formula with its assumptions documented).
