# RFC 0010: Launch-vehicle simulation track

- **Status:** Draft
- **Author(s):** Valenx contributors
- **Created:** 2026-05-29
- **Discussion PR:** _TBD_
- **Tracking issue:** _TBD_

---

## Summary

This RFC defines the full **launch-vehicle / spaceflight simulation
track** for Valenx and sequences it into discrete, independently-shippable
"bricks." It records what is **done today** (`valenx-astro`: native
atmosphere, gravity, drag, propulsion, ascent + closed-loop orbital
insertion), and lays out — honestly and concretely — every remaining
brick required to approach a *commercial-grade* rocket-design capability
(the "SpaceX could use it" bar), each with a scope, a definition of done,
and a rough effort estimate.

Its purpose is to make the path **legible**: so that "finish everything"
becomes a finite, ordered backlog rather than an open-ended aspiration,
and so each brick can be picked up, built, tested and merged on its own.

---

## Motivation

Valenx began as a generalist simulation shell (mostly external-solver
adapters, heavily bio/chem). It had **zero** rocket-specific physics. Two
bricks now exist (see below). Users keep asking, reasonably, "how close
is this to a tool SpaceX could design rockets with?"

The honest answer is: **very far** — a flight-grade suite is thousands of
engineer-years and a formal verification-and-validation (V&V) regime
against flight data. But "far" is not "formless." Every capability on the
way there is a known, bounded piece of engineering. This RFC writes them
all down and orders them, so progress is measurable and no one mistakes a
v1 brick for flight-certified software.

If we don't do this, the track grows ad hoc, scope is repeatedly
relitigated, and the gap to "commercial" stays vague.

---

## Guide-level explanation

The track is organised into **phases**. Each brick is a crate or module
with its own tests and a clear definition of done (DoD). A brick is only
"done" when it is implemented, unit-tested, documented, and (where
physics allows) validated against a known analytic result or published
reference.

### Phase 0 — Foundations *(DONE)*

| Brick | State | What it delivers |
|---|---|---|
| `valenx-astro` core | ✅ done | Planar 3-DOF point-mass dynamics, WGS-84 inverse-square gravity, US Standard Atmosphere 1976, Mach-dependent drag vs co-rotating air, RK4, Keplerian orbit determination. |
| Open-loop ascent | ✅ done | Vertical-rise → pitch-kick → gravity-turn guidance, multi-stage staging, max-Q / peak-g; reaches a bound orbit. |
| Closed-loop insertion | ✅ done | Ascent → coast-to-apoapsis → circularisation burn; reaches a near-circular ~300 km LEO (e ≈ 0.0003). |
| Native propulsion | ✅ done | Ideal-rocket (de Laval) engine model: c*, Γ, area-ratio→exit-Mach, thrust + Isp at any ambient pressure; `Stage::from_engine`. |
| 3-D orbital mechanics + J2 | ✅ done | `orbit3d`: classical orbital elements (rv↔COE), two-body + J2 RK4 propagator, secular J2 rates validated against the propagator. |
| Orbital maneuver planner | ✅ done | `maneuver`: Hohmann + bi-elliptic transfers and plane-change Δv budgets, validated against textbook values. |
| Lambert solver | ✅ done | `lambert`: universal-variable two-body targeting (rendezvous/intercept basis), validated by round-trip against the propagator. |
| Launch geometry | ✅ done | `launch`: azimuth↔inclination from latitude, minimum reachable inclination, Earth-rotation velocity bonus. |
| Ground tracks | ✅ done | `groundtrack`: sub-satellite lat/lon over the rotating Earth; max latitude = inclination, equatorial stays on equator. |
| Δv / propellant budgeting | ✅ done | `budget`: rocket-equation propellant ↔ Δv and mission-sequence feasibility / margin. |
| LLM/MCP control | ✅ done | `simulate_ascent`, `hohmann_transfer`, `launch_azimuth` MCP tools. |

### Phase 1 — Higher-fidelity flight mechanics

| Brick | Scope | DoD | Effort |
|---|---|---|---|
| 3-D orbital mechanics | 3-D ECI state, classical orbital elements (a, e, i, Ω, ω, ν) with exact rv↔COE conversion. | ✅ **done** (`orbit3d`): round-trip stable, inclination recovered, energy conserved over an orbit. |
| J2 gravity | Oblateness perturbation + secular nodal-regression / apsidal-precession rates. | ✅ **done** (`orbit3d::j2_accel`): propagated nodal regression matches the analytic rate to <5%; sun-sync inclination check passes. *(Higher-order EGM harmonics still open.)* |
| 3-D powered ascent | Embed the planar ascent into a 3-D orbit via launch azimuth → target inclination from a launch latitude. | ✅ **done** (`flight3d`): launch hits the target inclination exactly (`cos i = cos φ sin β`); the 3-D orbit shape matches the planar result. *(In-plane dynamics are planar; a native 3-D powered integrator with out-of-plane steering is the follow-up.)* |
| Rotating-atmosphere winds | Altitude-varying wind profiles; wind-relative drag. | ✅ **done** (`wind`): None / Constant / Gaussian-jet models raise max-Q and perturb the trajectory; zero-wind reproduces Phase 0 exactly. |
| 6-DOF rigid body — rotational core | Euler's equations + quaternion attitude under applied torque. | ✅ **done** (`rigidbody`): torque-free energy + inertial angular-momentum conservation and the analytic axisymmetric precession rate are all validated. |
| 6-DOF closed-loop demonstrator | Couple rotation to translation; thrust along the body axis; a PD pointing controller steers the attitude. | ✅ **done** (`flight6dof`): the loop slews to a commanded pointing, steers the thrust, and rejects a steady disturbance with bounded error — validated by control-theoretic behaviour. |
| 6-DOF production flight GNC | Guidance + navigation/estimation + control with actuator models (lag, saturation), aero moments, flex/slosh, validated against flight data. | A flight-representative vehicle flies a full validated mission. | XL — research-grade |

---

> ## ⛔ Boundary: tractable bricks vs. research-grade projects
>
> Everything marked ✅ above (Phase 0 + the analytic Phase-1 bricks) is
> **session-completable**: each has a closed-form or reference answer to
> validate against, so it can be built and proven correct in a single
> focused effort. **That track is now complete** — every analytic
> Phase-0 / Phase-1 brick is built and validated. The only Phase-1 item
> left, 6-DOF rigid-body flight, is research-grade (see below).
>
> **Phases 2–6 below are NOT session-completable.** Each is an
> independent, multi-month-to-multi-year engineering or research project
> (finite-rate combustion, compressible/hypersonic CFD, 6-DOF GNC,
> aeroheating/TPS, coupled loads, and a formal validation-against-
> flight-data regime). They cannot be honestly "finished" by an
> automated loop: there is no ground-truth oracle to validate generated
> code against, so producing code here without a dedicated project and
> real validation data would yield plausible-looking but unverified —
> and therefore unsafe — results. They are documented here as the
> roadmap, to be picked up as funded, staffed work items, **not** as
> tasks to be auto-generated.

### Phase 2 — Propulsion depth *(research-grade — see boundary note above)*

| Brick | Scope | DoD | Effort |
|---|---|---|---|
| Finite-rate combustion | Replace frozen-gas assumption with equilibrium/finite-rate chemistry (drive `valenx-qchem`/Cantera). | Chamber c* / Isp match CEA/RPA for a reference propellant within a few %. | L |
| Engine cycles + feed system | Gas-generator / staged-combustion / electric-pump cycles, tank pressures, turbopump power balance, throttling, ullage. | Predicts throttle-down Isp and feed-pressure limits for a reference engine. | L |
| Nozzle separation + altitude comp | Summerfield separation criterion, dual-bell / aerospike. | Flags flow separation regimes; aerospike altitude-compensation curve is monotone-correct. | M |

### Phase 3 — Aerothermodynamics & reentry

| Brick | Scope | DoD | Effort |
|---|---|---|---|
| Compressible CFD | Native compressible/transonic→supersonic solver, or a validated SU2/OpenFOAM-compressible pipeline + body-fitted meshing. | Recovers a known wedge/cone shock angle and Cd vs Mach for a standard body. | XL |
| Hypersonic + real-gas | High-temperature real-gas effects, equilibrium/non-equilibrium chemistry behind the shock. | Stagnation heating matches Fay–Riddell for a sphere within engineering tolerance. | XL |
| Aeroheating + TPS | Surface heat flux, ablation, thermal-protection sizing. | TPS recession for a reference reentry matches a published case. | L |

### Phase 4 — Guidance, navigation & control

| Brick | Scope | DoD | Effort |
|---|---|---|---|
| Powered explicit guidance (PEG) | Closed-loop optimal ascent targeting (replaces the open-loop kick-tuning). | Hits a commanded orbit (a/e/i) across a range of vehicles without per-case tuning. | L |
| Trajectory optimisation | Pseudospectral / collocation optimal-control for ascent, RTLS/ASDS booster return, landing burns. | Reproduces a published minimum-propellant ascent or a hoverslam landing. | XL |
| Navigation + estimation | IMU/GNSS sensor models, Kalman filtering, dispersions, Monte-Carlo. | 1000-run Monte-Carlo produces a sane insertion-accuracy ellipse. | L |

### Phase 5 — Coupled multiphysics & loads

| Brick | Scope | DoD | Effort |
|---|---|---|---|
| Aero-thermal-structural coupling | Couple `valenx-aero`/CFD ↔ `valenx-fem` ↔ thermal via the existing preCICE path. | A flight-loads case converges and matches a benchmark. | XL |
| Structural loads & margins | Max-Q gust loads, stage-separation, pogo, slosh, flutter, factor-of-safety reporting. | Produces a loads envelope and margins for a reference airframe. | L |

### Phase 6 — Productisation & V&V

| Brick | Scope | DoD | Effort |
|---|---|---|---|
| Formal V&V regime | Documented verification (method-of-manufactured-solutions, convergence studies) + validation against flight/test data. | A V&V report exists per solver with quantified error bars. | XL |
| Mission/UX front end | Vehicle builder, mission planner, trajectory + telemetry visualisation in the desktop shell, end-to-end LLM-driven workflows. | A non-expert can build a vehicle and fly a mission from the GUI or by chat. | L |
| QA, traceability, support | Requirements traceability, certification-grade testing, long-term support posture. | Meets the process bar a flight programme requires. | XL |

Effort key: **S** ≈ days, **M** ≈ weeks, **L** ≈ months, **XL** ≈ many
months to years (often a sub-project in its own right).

---

## Reference-level explanation

- New physics lands as new crates (`valenx-astro` siblings) or modules,
  never by weakening an existing solver's documented scope.
- Every brick keeps the project's **honest-scope** convention: the crate
  docs state precisely what is and is not modelled, and what each
  omission would take to close (mirroring `valenx-aero` and
  `valenx-astro`).
- LLM/MCP exposure is added per brick as its native crate ships, the same
  way `simulate_ascent` was added.
- No brick is marked "done" without tests; physics bricks additionally
  require a validation check against an analytic result or published
  reference.

---

## Drawbacks

- The full track is enormous; most bricks are independently the size of a
  serious open-source project. Marketing this as "near commercial" before
  the V&V phase would be dishonest and is explicitly out of bounds.

## Rationale and alternatives

- **Why phase it:** each brick is independently useful (e.g. 3-D + J2 +
  PEG already makes a credible mission-design tool long before
  hypersonics or 6-DOF).
- **Alternative — wrap existing tools only:** adapters to GMAT / OpenRocket
  / RPA would get capability faster but reintroduce the external-tool
  dependency the native engines deliberately avoid. The two are
  complementary; adapters can fill gaps while native bricks mature.

## Prior art

NASA GMAT, Astrogator, POST2, OpenRocket, RocketPy, RPA, NASA CEA,
SU2 (compressible CFD), preCICE (coupling). Each maps to one or more
bricks above.

## Unresolved questions

- Native compressible CFD vs. a hardened SU2 pipeline for Phase 3.
- How far to push 6-DOF before it earns its (XL) cost.
- What validation datasets are publicly available for the V&V phase.

## Future possibilities

Interplanetary trajectory design (patched conics → full n-body),
launch-site / range-safety modelling, constellation/coverage analysis,
reusable-booster economics.
