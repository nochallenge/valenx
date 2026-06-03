# RFC 0011: Neural-interface / BCI simulation track (`valenx-neuro`)

- **Status:** Implemented (v1) — see Amendments
- **Author(s):** nochallenge
- **Created:** 2026-06-02
- **Discussion PR:** _TBD_
- **Tracking issue:** _TBD_

---

## Summary

Introduce `valenx-neuro`, a native-Rust neural-interface simulation suite
modelling the physics of an implanted stimulating electrode and its effect on
nearby neurons — the core physics behind a brain–computer interface (BCI). Five
coupled modules:

1. **Extracellular field** — the electric potential an electrode sets up in
   tissue (FEM).
2. **Hodgkin–Huxley cable** — the membrane dynamics of axons.
3. **Activating function** — how the field drives the membrane (the coupling).
4. **Bioheat** — tissue temperature rise from stimulation (Pennes).
5. **Electrode–tissue impedance** — what the electrode "sees."

Each module is validated against a closed-form or textbook result, the same way
`valenx-astro` is checked against the Hohmann and J2 results. A `valenx-app`
workbench provides setup, background compute, 3-D playback, and plots, mirroring
the reaction-dynamics workbench.

---

## Motivation

Valenx spans aerospace, engineering, chemistry, and biology, but has **no
neural-engineering physics**. Neural-interface modelling — the domain of
cochlear implants, deep-brain stimulation, and high-density electrode arrays — is
a well-defined, physically grounded simulation problem: given an electrode, a
stimulation waveform, and tissue, *which neurons fire, how much charge and heat
is delivered, and what does the electrode measure?*

These questions have established models with analytic checks — the **Rattay
activating function**, the **Hodgkin–Huxley** equations, the **Pennes bioheat**
equation — which makes them a clean fit for Valenx's "native solver, validated
against a reference" pattern. The work also extends the suite into
neuroengineering, a deliberate target.

**Honest scope.** This is a **research / education-grade** first-principles
model, not clinical or device-certification software. It uses quasi-static
fields, idealized (procedural) tissue and point/disk electrode geometry, and
standard membrane models. It is explicitly **not** a validated neurostimulation
design tool and must not be represented as one.

---

## Guide-level explanation

The canonical workflow:

1. Place an **electrode** in a block of **tissue** (conductivity σ).
2. Define a **stimulation pulse** — amplitude (µA), pulse width (ms),
   cathodic / anodic / biphasic.
3. Lay down one or more **axons** nearby (position, diameter).
4. Press **Run**.

Valenx then solves the extracellular potential field, samples it along each
axon, computes the activating drive, integrates the Hodgkin–Huxley membrane
dynamics over time, and reports:

- **which axons fired**, and the action potential propagating along them;
- the **recruitment curve** — fraction of axons recruited vs. stimulus amplitude;
- the **tissue temperature rise** from the deposited power;
- the **electrode impedance spectrum**.

The 3-D viewport animates the action potential travelling along recruited axons
with the extracellular field shown as a heatmap on a slice plane.

The intuition this teaches is the heart of BCI stimulation design: a **cathodic**
pulse is more efficient at firing a nearby fiber than an **anodic** one, and the
current needed to recruit a fiber grows with the **square of its distance** from
the electrode. The suite reproduces both from first principles.

---

## Reference-level explanation

### Crate layout

A new top-level crate `valenx-neuro` (mirrors `valenx-reactdyn`):

```
crates/valenx-neuro/src/
  lib.rs        — public API + the Simulation orchestrator
  units.rs      — unit conventions + conversions (the one place numbers convert)
  error.rs      — typed errors
  scene.rs      — procedural geometry: tissue block, electrode, axon paths
  field.rs      — extracellular field: −∇·(σ∇φ)=I via valenx-fem, relabeled
  cable.rs      — Hodgkin–Huxley multi-compartment axon
  activating.rs — Rattay activating function: 3-D φ sampled along a fiber → drive
  bioheat.rs    — Pennes bioheat: valenx-fem conduction + perfusion + Joule source
  impedance.rs  — electrode model: access resistance + CPE double layer → Z(ω)
  engine.rs     — orchestrates a run → Trajectory (V(x,t), spikes, recruitment)
```

The `valenx-app` workbench lives in `crates/valenx-app/src/neuro_workbench.rs`,
wired into the View menu and draw dispatch exactly like `reactdyn_workbench.rs`.

### Units discipline

Unit mismatches are *the* classic bug in coupled electrophysiology (the
analogue of the AIMD unit bug). One convention, converted only in `units.rs`,
with a round-trip test:

| Quantity | Unit |
|---|---|
| Potential | mV |
| Time | ms |
| Current | µA |
| Conductivity σ | S/m |
| Tissue length | mm |
| Compartment length / fiber diameter | µm |
| Membrane capacitance C_m | µF/cm² |
| Conductances g | mS/cm² |
| Temperature | °C (ΔT in K) |

Every cross-module hand-off (field mV → cable extracellular drive) passes
through a typed conversion.

### Module 1 — Extracellular field (`field.rs`)

Quasi-static current conduction in tissue:

```
−∇·(σ ∇φ) = I_v
```

This is **identical in form** to the steady heat equation `−∇·(k∇T)=q` already
solved by `valenx-fem`'s `thermal_solver` (σ↔k, φ↔T, injected current↔heat
source). Module 1 reuses that tetrahedral assembly + sparse Cholesky factorization
and relabels units — a real FEM field, not a toy. The electrode is a
current-injection (Neumann) boundary; the far field is grounded (Dirichlet φ=0).

Because the problem is linear, the field is solved once for a unit current and
scaled by the waveform.

**Validation:** a point current source `I` in an unbounded homogeneous medium
gives `φ(r) = I / (4πσr)`. The FEM solution matches this away from the source
singularity, within mesh-discretization error.

### Module 2 — Hodgkin–Huxley cable (`cable.rs`)

A multi-compartment axon. Per compartment:

```
C_m dV/dt = −(I_Na + I_K + I_L) + I_axial + I_drive
I_Na = g_Na m³ h (V − E_Na)
I_K  = g_K  n⁴  (V − E_K)
I_L  = g_L      (V − E_L)
```

with the HH (1952) gating ODEs `dm/dt = α_m(V)(1−m) − β_m(V) m` (and `h`, `n`),
the standard rate functions, and squid-axon parameters (`g_Na=120`, `g_K=36`,
`g_L=0.3` mS/cm²; `C_m=1` µF/cm²; rest ≈ −65 mV). Axial coupling
`I_axial = (a / 2R_i) ∂²V/∂x²` ties compartments into a cable. Integrated with
fixed-step RK4 (small Δt; HH is mildly stiff — see Unresolved Questions).

**Validation:** a single compartment given a supra-threshold stimulus fires a
**textbook action potential** — ~100 mV overshoot from rest, a firing threshold,
and a refractory period. A multi-compartment cable produces a **propagating**
action potential with a **conduction velocity** in the expected range for the
chosen diameter.

### Module 3 — Activating function (`activating.rs`)

The coupling. The extracellular potential `V_e(x)` is sampled from Module 1's
field along each axon's centerline. With the membrane seeing `V_m = V_i − V_e`,
the cable equation gains a driving term — the **Rattay activating function** —
proportional to the *second spatial derivative* of the extracellular potential
along the fiber:

```
f_n ∝ ∂²V_e/∂x²   (per compartment, the discrete second difference)
```

**Validation:**
- **Sign / polarity:** under a **cathodic** electrode (I<0), `V_e<0` and
  `∂²V_e/∂x²>0` at the nearest node → **depolarizing** (excitatory); anodic is
  hyperpolarizing under the electrode with excitatory flanks. The **cathodic
  threshold is lower than the anodic** threshold.
- **Strength–distance:** the threshold current to fire a fiber scales as **r²**
  with electrode-to-fiber distance.

### Module 4 — Bioheat (`bioheat.rs`)

The Pennes bioheat equation, steady state first:

```
∇·(k ∇T) − ω_b ρ_b c_b (T − T_a) + Q = 0
```

Reuses `valenx-fem` conduction (same solver as Module 1), adds **perfusion** as a
reaction (diagonal) term, and a heat source `Q = σ|∇φ|²` — the resistive (Joule)
heating computed from Module 1's field.

**Validation:** a steady point heat source without perfusion gives
`ΔT = Q/(4πk r)` (the thermal twin of the electric point source); with perfusion,
`ΔT = Q/(4πk r)·e^(−r/L)` with penetration depth `L=√(k/ω_b ρ_b c_b)`. The FEM
result matches.

### Module 5 — Electrode–tissue impedance (`impedance.rs`)

A lumped electrode model. For a disk electrode of radius `a` on tissue of
conductivity σ, the **access (spreading) resistance** is:

```
R_a = 1 / (4 σ a)
```

In series with a **double-layer constant-phase element** `Z_CPE = 1/(Q(jω)^n)`,
optionally with a Faradaic charge-transfer resistance `R_ct` in parallel:

```
Z(ω) = R_a + ( R_ct ∥ Z_CPE )
```

**Validation:** `R_a = 1/(4σa)` against the analytic spreading-resistance
formula; the `|Z|(ω)` Bode shape (capacitive at low frequency, resistive plateau
at high frequency).

### Data flow (one run)

```
Scene (tissue σ, electrode pos/radius/waveform, axon bundle)
  → mesh tissue (structured tet grid v1)
  → field solve (unit current, linear → scale by waveform)        [Module 1]
  → sample V_e along each axon centerline                          [Module 3]
  → HH integration over time with the activating drive            [Modules 2+3]
  → record V(x,t), spike times, recruited fraction
  → sweep amplitude → recruitment curve
  → Joule power → steady ΔT field                                 [Module 4]
  → electrode geometry + σ → Z(ω) + voltage compliance            [Module 5]
  → Trajectory + summaries → workbench
```

### Workbench (`neuro_workbench.rs`)

View-menu toggle. Setup panel (tissue preset — gray/white matter σ; electrode
waveform — amplitude / pulse width / mono- or biphasic; axon count / diameter /
layout) → **Run** in a background thread → 3-D viewport (tissue box, electrode
glyph, axons colored by membrane potential, field heatmap on a slice plane, a
time scrubber animating the AP) → plots (membrane V(t), recruitment curve, ΔT,
`|Z|(ω)` Bode). Reuses the reactdyn 3-D playback scaffolding and `egui_plot`.

---

## Validation targets

| Module | Reference | Check |
|---|---|---|
| Extracellular field | point source, homogeneous medium | φ = I/(4πσr) within mesh error |
| HH cable | Hodgkin–Huxley 1952 | AP ~100 mV overshoot, threshold, refractory period; conduction velocity in range |
| Activating function | Rattay 1986 | cathodic threshold < anodic; recruitment threshold rises with electrode distance |
| Bioheat | Pennes, point source | ΔT = Q/(4πk r)·e^(−r/L) vs analytic |
| Impedance | disk electrode | R_a = 1/(4σa); capacitive→resistive Bode shape |

---

## Build phases (definition of done)

Strict TDD order — each green and committed before the next:

0. `units.rs` round-trip tests.
1. HH **single compartment** → textbook AP (no field).
2. HH **multi-compartment cable** → propagating AP + conduction velocity.
3. **Extracellular field** (reuse `valenx-fem`) → point-source φ.
4. **Activating function** → cathodic/anodic sign + strength–distance r².
5. **Coupled stimulation** → electrode current recruits the nearest axon at the
   expected threshold; recruitment curve.
6. **Bioheat** → analytic point-source ΔT.
7. **Impedance** → R_a = 1/(4σa) + Bode shape.
8. **Workbench** wiring → headless UI tests (panel draws, Run path, bad-input
   handling), mirroring the reactdyn workbench tests.

A brick is "done" when implemented, unit-tested, documented, and validated
against its reference.

---

## Drawbacks

- **Research-grade, not clinical.** Idealized geometry and standard membrane
  models; not a substitute for validated neurostimulation design software, and
  must never be presented as one.
- **Quasi-static fields.** No electromagnetic wave propagation — valid below
  ~100 kHz, which covers neural stimulation, but it is an assumption.
- **Idealized geometry in v1.** Homogeneous, procedurally meshed tissue; point /
  disk electrodes. Heterogeneous/anisotropic tissue and realistic electrode
  arrays are later bricks.
- **Compute-then-visualize.** The coupled FEM + cable solve is seconds-to-minutes
  for a bundle, not interactive/real-time.
- **Squid HH first.** Mammalian / myelinated fiber models (CRRSS, MRG) with
  saltatory conduction are a later brick.

---

## Alternatives considered

- **Dedicated EM solver inside `valenx-neuro`.** Rejected: duplicates the
  validated `valenx-fem` solver for more work and more risk. The PDE is identical;
  reuse wins.
- **Analytic-only field** (`φ = I/4πσr`, no mesh). Rejected: it is not the *FEM
  extracellular field* the suite advertises and cannot represent electrode
  geometry, tissue boundaries, or heterogeneity. Kept instead as the *validation
  reference* for the FEM.
- **Full 3-D bidomain.** Rejected: overkill for extracellular stimulation of
  passive-environment fibers; large and slow with no benefit at this scope.
- **Real-time interactive solve.** Rejected: the coupled FEM + HH integration is
  not interactive-cheap; reactdyn-style compute-then-playback is the right fit.

---

## Unresolved questions

- **Mesh generation:** reuse `valenx-fem`'s `meshgen`, or a simple structured
  tetrahedral grid for v1? *Lean:* structured grid for v1, swap later.
- **Cable integrator:** fixed-step RK4 vs. adaptive vs. implicit (HH is mildly
  stiff). *Lean:* RK4 with a conservative Δt for v1; revisit if stability bites.
- **Myelinated fibers:** v1 is unmyelinated HH only; node-of-Ranvier / saltatory
  models deferred to a later brick.
- **Bioheat / impedance depth in v1:** steady-state bioheat + R+CPE impedance
  first; transient bioheat and full EIS spectra later.
- **Biphasic charge-balancing** waveform details (inter-phase gap, asymmetry).

---

## Amendments

**2026-06-02 — v1 implemented (`valenx-neuro`, 20 tests).** Built per this RFC:
extracellular FEM field (point-source φ = I/4πσr validated), Hodgkin–Huxley
compartment + cable (textbook action potential + propagation), the Rattay
activating function (cathodic/anodic sign), coupled recruitment, bioheat, and
electrode impedance, plus a `valenx-app` workbench. Honest research-grade
scoping notes / deviations from the original targets:

- **Bioheat is conduction-only.** `valenx-fem`'s solver has no reaction term, so
  v1 solves `−∇·(k∇T)=Q` (point-source ΔT validated); the Pennes perfusion term
  (`e^{−r/L}` cooling) is deferred.
- **1-mm cable compartments.** Explicit RK4 axial integration is stable only for
  roughly `αΔt ≤ 0.5`; at 100 µm it diverges in the sub-threshold (linear)
  regime. v1 uses 1-mm compartments (stable; the ~6 mm space constant is well
  resolved). A finer cable needs an implicit (Crank–Nicolson) axial solver — the
  proper fix.
- **Strength–distance rises with distance, not exactly ∝ r².** On the coarse
  (2-mm) FEM field the measured threshold exponent over 0.5–2 mm is ≈ 0.7–1, not
  the textbook current-distance `r²`; the qualitative law (deeper fibers need
  more current) holds, and a finer field should steepen it.

**2026-06-03 — Phase 2 implemented (research-grade upgrades; `valenx-neuro` now
33 tests).** Five modules added on the v1 base; the v1 scoping notes above are
largely **resolved**:

- **Implicit cable solver** (`membrane`) — backward-Euler axial diffusion with a
  Thomas tridiagonal solve, behind a generic `Membrane` trait. Unconditionally
  stable: the exact 100 µm sub-threshold case that diverged the v1 explicit RK4
  to +∞ now stays bounded and still propagates. **Resolves the "1-mm cable
  compartments" note** above.
- **Myelinated mammalian fiber** (`myelinated`) — active Hodgkin–Huxley nodes of
  Ranvier joined by near-transparent myelinated internodes (length ∝ diameter).
  Conduction velocity reproduces the empirical **CV ≈ 6·D** rule within ~6%
  (57 / 113 m/s at 10 / 20 µm) and scales ∝ D (saltatory), not ∝ √D.
- **Strength–duration** (`strength_duration`) — rheobase + chronaxie (1.65 ms,
  ≈ ½ the membrane time constant) by bisection; the Lapicque/Weiss
  constant-charge law holds to < 1% at short pulse widths.
- **Anisotropic / heterogeneous FEM field** (`aniso_field`) — a from-scratch
  solve of `−∇·(σ∇φ)=I` with a per-element 3×3 conductivity tensor, by conjugate
  gradient. Validated against the closed-form anisotropic point source
  `φ = I/(4π√(detσ)·√(rᵀσ⁻¹r))` to within ~10% (under 4% away from the source).
- **Multi-contact steering** (`steering`) — field superposition (validated to
  solver tolerance) and electronic current steering that shifts the stimulation
  focus without moving the lead.

Honest ceiling: this is **real neurostimulation-research-grade** modelling
(standard membrane models, idealized geometry, quasi-static fields), not
Neuralink-production or clinical software — that needs their hardware, data, and
regulatory pipeline, not just code.

**2026-06-03 — Extracellular recording added (forward-EAP, `recording`).** The
read-out side of a neural interface (the rest of the crate stimulates). A firing
axon's transmembrane currents act as extracellular point sources,
`φ_e = 1/(4πσ_e)·Σ_k I_m,k/r_k`; the membrane currents are built from the axial-
current divergence, so they conserve charge exactly (`|Σ I_m| ≈ 1e-21 A`).
Validated: the recorded waveform is **biphasic with a dominant-negative (sink)
phase** — the textbook extracellular action potential — and its amplitude falls
off with electrode distance. `valenx-neuro` now 36 tests.
