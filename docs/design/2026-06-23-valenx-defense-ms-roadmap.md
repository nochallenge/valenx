# valenx Defense / Modeling-&-Simulation (M&S) Capability Roadmap

**Date:** 2026-06-23
**Status:** queued for the autonomous build loop (branch `feat/capabilities`), after / alongside the COLMAP-SfM track.

---

## 0. Scope & honesty boundary (READ FIRST)

valenx will build **research-grade, in-house, AI-drivable analogs of the *public* capability
classes** of defense modeling & simulation — the same **dual-use engineering posture** as Ansys,
MATLAB, MSC, Siemens, STK, and the DoD CREATE design suites. valenx is a *defense-adjacent
engineering & simulation platform*, **not** a weapons program.

This boundary is a hard gate on every track below:

- **General engineering & simulation — built fully.** CFD, FEM, flight dynamics, orbital mechanics,
  uncertainty quantification, digital engineering, sensor/EM physics. These are dual-use; defense is
  one customer alongside civilian aerospace, space, automotive, and telecom.
- **Dual-use physics — framed neutrally.** Survivability/protection (not lethality), vehicle GNC
  (not weapon targeting), signature *design/analysis* (not malicious detection-evasion). We build the
  neutral physics and the *defensive / design* application.
- **Excluded.** Dedicated weapons-lethality, warhead/penetration "kill" modeling, fire-control /
  targeting, and anything whose primary purpose is to harm more effectively. Where a physics domain
  is shared (impact, blast, EM), we build the **defensive / design** side only — armor & hardening,
  link budgets, low-observable *design* — never offensive optimization.
- **Research-grade, validation-pending.** Honest engineering increments, **not** accredited (VV&A)
  defense systems. Where the real-world gap is *validated* tooling (e.g. hypersonics), we build the
  solver and label validation as **pending** — the ground-truth data is not in hand and we will not
  fake it.

Every track also obeys the standing valenx gates:

1. **Reactive front-end + AI-drivable** — each capability lands a workbench *and* is agent-openable
   (registered under a workbench id the agent bridge can drive, named widgets). Not a headless crate.
2. **TDD + gated green** — `build` · scoped `test` · `clippy --all-targets -D warnings` · `fmt` ·
   `rustdoc -D warnings`.
3. **Code + security review to zero findings** before advancing.
4. **In-house** — port permissive (BSD/MIT/Apache/MPL) directly into Rust with attribution in
   `THIRD-PARTY-NOTICES`; never port GPL/copyleft (clean-room from published math instead).
5. **BLOCKED, not faked** — anything that needs data/weights/licensed inputs we don't have is
   flagged honestly, never stubbed as if working.

All numbers research-grade unless cross-validated against named ground truth (see `docs/VALIDATION.md`).

---

## 1. The queue — nine tracks

Each track is mapped to the user's three buckets — what the military **HAS** (mature tools), what it
is **BUILDING** (active investment), and the real **GAP** (where the pull is) — then to what valenx
can honestly build and the dual-use posture.

### M1 — UAS / counter-UAS design → sim → trade study
- **Has:** commercial drone CAD; hobby/industry flight sims.
- **Building:** Collaborative Combat Aircraft / "loyal wingman", drone swarms, counter-UAS.
- **Gap:** a *fast, iterative* small-UAS design-and-simulate loop with trade studies.
- **valenx leverage (~70% already in-house):** `valenx-rotor` BEMT (just shipped), `flight6dof`,
  drag/aero, battery, propulsion, CFD.
- **Build:** multirotor & fixed-wing **vehicle assembly** → integrated performance (thrust, hover
  endurance, range, payload) → **mission/trade studies** (sweep design params, pareto fronts).
  Counter-UAS = detection/intercept **geometry & timeline** analysis (defensive).
- **Dual-use posture:** civilian drone design is the same tool; counter-UAS framed as defensive
  detect/track/intercept *geometry*, not weapon employment.

### M2 — Mission / engagement constructive-simulation framework (AFSIM-class)
- **Has:** AFSIM, OneSAF, JCATS (government-owned engagement & wargame sims).
- **Building:** JADC2 / multi-domain, kill-web analysis, autonomy experimentation.
- **Gap:** interoperable, AI-drivable, fast mission-level analysis.
- **Build:** a **discrete-event / agent simulation framework** — entities, movers (reuse
  `flight6dof` / `valenx-astro`), sensors (detection-range models), comms/networks, scheduler.
  Engagement outcomes stay **abstract & probabilistic** (probability-of-kill as an *input parameter*,
  Lanchester-style attrition) — **no** detailed lethality or targeting.
- **Dual-use posture:** the framework is general (also logistics, epidemiology, traffic, policy
  wargaming used by think-tanks/universities). We build infrastructure + analysis, not kill chains.

### M3 — Digital-engineering / MBSE / digital-twin spine
- **Has:** commercial MBSE (Cameo/Magic, DOORS).
- **Building:** the DoD Digital Engineering Strategy — digital thread, model-based acquisition,
  sustainment digital twins.
- **Gap:** a *real* digital thread tying coupled multiphysics into one design environment.
- **Build:** extend `valenx-orchestrator` into a **systems-model + requirements-trace + trade-study +
  digital-thread** layer that wires valenx's physics solvers into coherent, versioned design studies
  with traceability. Pure systems-engineering software; zero dual-use concern.

### M4 — Uncertainty quantification + surrogates + sensitivity  *(cross-cutting enabler — build first)*
- **Gap:** decision-relevant UQ *everywhere* — leaders need confidence bounds, not point estimates.
- **Build:** a `valenx-uq` crate wrapping **any** valenx solver: Monte-Carlo & Latin-hypercube
  sampling, polynomial-chaos expansion, Sobol/Morris global sensitivity, surrogate models
  (Gaussian-process / polynomial / small NN), model calibration, and propagated confidence intervals.
- **Dual-use posture:** general numerical methods; the single highest-leverage non-weaponized item,
  and a force-multiplier for every other track.

### M5 — Hypersonic / high-enthalpy aerothermodynamics  *(research-grade; validation is the honest gap)*
- **Has:** CREATE-AV, commercial CFD.
- **Building:** hypersonic vehicle design — aerothermal heating, GNC, thermal-protection materials.
- **Gap:** **validated** hypersonic / high-enthalpy CFD (ground-test capacity can't cover the regime).
- **Build:** extend valenx CFD toward high-Mach compressible flow, aerodynamic heating, and
  ablation / TPS (thermal-protection-system) sizing. **Honest frame:** vehicle design (also civilian
  reentry & space access); research-grade, **validation pending** (no ground-truth data — labeled, not faked).

### M6 — Space domain / orbital operations & SSA (STK-class)
- **Has:** STK (ubiquitous).
- **Building / Gap:** AI-drivable space mission design; space situational awareness at scale.
- **Build:** extend `valenx-astro` — access/coverage analysis, constellation design, conjunction
  screening & debris (SSA), maneuver planning, ground-station scheduling.
- **Dual-use posture:** civilian + defense space; the analysis is the same either way.

### M7 — Survivability & protection (defensive blast / impact)
- **Has:** LS-DYNA, CTH (Sandia), ALE3D (LLNL).
- **Gap:** fast survivability trade-offs.
- **Build:** structural response to blast & impact **for protection** — armor sizing, hardened
  structures, vehicle/occupant survivability — using valenx FEM / explicit dynamics.
- **Dual-use posture:** **defensive only.** Same physics as civil crash safety and blast-resistant
  building design. We build "how to *protect*", never "how to *penetrate*". No warhead/lethality models.

### M8 — Sensor / RF / EM & signature *design*
- **Has:** Xpatch (RCS), commercial EM.
- **Building / Gap:** integrated, fast signature-in-the-loop design.
- **Build:** EM propagation, radar detection-range / link budgets, and radar-cross-section prediction
  (physical optics) for **design & analysis** — antenna design, low-observable *design*, sensor
  trade studies.
- **Dual-use posture:** shared with telecom & radar engineering; framed as design/analysis, never as
  malicious detection-evasion.

### M9 — Autonomy V&V / T&E methodology
- **Gap:** trusted verification & validation of learned / autonomous systems — the #1 chokepoint for
  fielding autonomy.
- **Build:** a test-harness + scenario-generation + coverage-metrics + assurance-evidence framework
  for autonomous / AI components (also general AI safety). Aligns directly with valenx's
  AI-drivable-first thesis.
- **Dual-use posture:** methodology & tooling; no dual-use concern.

---

## 2. Priority order

Lead with the cross-cutting enabler and valenx's existing strengths, then breadth:

1. **M4** — UQ + surrogates (enables every other track; zero dual-use concern).
2. **M1** — UAS / counter-UAS (valenx already owns ~70%; real, current gap).
3. **M3** — digital-engineering spine (extends the orchestrator; ties it together).
4. **M6** — space / SSA (extends `valenx-astro`).
5. **M2** — mission-sim framework.
6. **M5** — hypersonic aerothermo (research-grade).
7. **M7** — survivability.
8. **M8** — sensor / RF / signature.
9. **M9** — autonomy V&V.

Order may interleave where crates are independent and parallelizable.

---

## 3. Mechanism

Same autonomous loop, same session, spawning subagents:

> build subagent (TDD + gated green) → orchestrator commits email-safe (`nochallenge
> <201502404+nochallenge@users.noreply.github.com>`, leak-scan 0, stage only real files) →
> code + security review to **zero findings** → reactive + AI-drivable workbench → next.

Runs after / alongside the COLMAP-SfM track. **GitHub held (local commits only)** until an explicit
"push". Each track is its own crate(s) + workbench so reviews stay scoped.
