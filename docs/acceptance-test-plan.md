# valenx User-Acceptance Test Plan (AI-driven walkthrough)

**Goal:** exercise valenx end-to-end *as a user would* — open a workbench, set inputs, run,
check the result — across the whole breadth of the app. Every test is expressed as a sequence
the **agent file-bridge** can drive (so it doubles as a proof that the app is AI-drivable), and
each is equally runnable by hand in the GUI.

## How we run these
- **Drive:** write the listed `{cmd:...}` JSONL lines to the bridge inbox; valenx polls (~1 Hz)
  and applies them. `OpenWorkbench{id}` opens a panel · `SetControl{name,value}` sets a labelled
  param by its on-screen caption · `RunCommand{id}` fires a palette action · `set_view/add_sketch_point/...`
  drive the viewport. `ListControls{workbench}` prints the exact settable captions for discovery.
- **Prereq:** a valenx build from `feat/capabilities` (the bridge commands are new) with
  `VALENX_ASSISTANT_INBOX/FEED` set. Relaunch needed before the AI-driven run.
- **Verify:** inputs + run are fully bridge-driven. Reading the *computed number* back is currently
  visual (the panel readout) — **likely our first gap to close** is a `ReadReadout`/result-to-feed
  bridge command so pass/fail can be auto-checked. Until then the human reads the panel.
- **Failure protocol (per the user):** if a step can't be driven, or the result is wrong/missing/crashes
  → log it → add the code to valenx → re-run → continue. Failures are the point.

Legend: **Goal · Steps · Expect · Pass**. `set "<caption>" = <v>` = `SetControl`.

---

## A. Aerospace & flight
- **A1 Rocket ascent.** Open `rocket` → set propellant/mass/thrust params → `RunCommand run.case`.
  Expect: a 3-D launch vehicle + an ascent trajectory/altitude plot. Pass: trajectory renders, apogee finite.
- **A2 UAS endurance.** Open `uas` → set mass, disk area, battery Wh, payload → run.
  Expect: hover power, endurance, range, payload margin. Pass: endurance ≈ E/P (minutes), no panic.
- **A3 Counter-UAS intercept.** Open `uas` → set threat track + interceptor speed + sensor range.
  Expect: time-to-intercept + intercept-point plan view. Pass: head-on TTI ≈ range/closing-speed; faster-opening target → "no intercept".
- **A4 Aero wind tunnel.** Open `aero` → set test body, AoA, speed → run.
  Expect: CL/CD + a flow-field view. Pass: CL rises with AoA below stall; CD finite.
- **A5 Engine cycle.** Open `engine` → set chamber pressure / mixture → run.
  Expect: combustion temps + cycle power balance. Pass: closes to a sane chamber pressure, no NaN.

## B. Astro & space
- **B1 Hohmann transfer.** Open `astro` → planner → set from/to altitude → run. Expect: Δv1, Δv2, transfer time. Pass: Δv ≈ vis-viva values.
- **B2 Reentry heating.** Open `astro` → reentry → set entry velocity/angle. Expect: peak deceleration + peak heat flux + altitude. Pass: matches Allen–Eggers/Sutton–Graves ballpark.

## C. CAD & geometry
- **C1 Sketch→extrude.** Open `cad` → `add_sketch_point` ×4 (a square) → `extrude_sketch{height:10}`.
  Expect: a 3-D box part appears. Pass: solid renders; volume ≈ area×height.
- **C2 2-D draft.** Open `draft2d` → `add_2d_line` + `add_2d_circle`. Expect: entities drawn. Pass: entity count increments, both visible.
- **C3 Sheet-metal bend.** Open `sheetmetal` → set thickness, bend angle, radius. Expect: flat pattern + bend allowance. Pass: allowance > 0, scales with angle.
- **C4 Mesh transform.** Open `meshtoolbox` → set translate/scale/rotate → apply. Expect: mesh moves. Pass: AABB shifts by the set amount.
- **C5 Surface revolve.** Open `cad`/surface → NURBS surface of revolution. Expect: a swept surface. Pass: renders, closed where expected.

## D. Structural & FEM
- **D1 Beam FEM.** Open `fem` → set load + geometry → run static + modal. Expect: deflection + stress + mode frequencies. Pass: tip deflection ≈ PL³/3EI; f1 > 0.
- **D2 Survivability.** Open `survivability` → set charge W, standoff R, plate → run.
  Expect: Friedlander P(t), P-I diagram, armor thickness, occupant margin. Pass: overpressure ↓ with R; armor ↑ with threat KE.
- **D3 Frames/truss.** Open `frames` → set member sizes/loads → run. Expect: member forces + deflection. Pass: equilibrium holds (ΣF≈0).
- **D4 Pressure vessel.** Open `pressurevessel` → set P, R, t. Expect: hoop/longitudinal stress. Pass: σ_hoop ≈ PR/t.
- **D5 Reinforcement.** Open `reinforcement` → set beam + rebar → run. Expect: 3-D rebar layout + capacity. Pass: renders, Mn finite.

## E. Fluids & CFD
- **E1 CFD cavity.** Open `cfd` → set lid speed, Reynolds, grid → run. Expect: velocity field + residual plot. Pass: residuals fall; recirculation visible.
- **E2 SPH.** Open `fluids` → set particle count, dt, steps → run. Expect: an SPH free-surface flow. Pass: particles stay bounded, density floored (no blow-up).
- **E3 Ocean.** Open `ocean` → set wave amplitude/period + hull → run. Expect: Gerstner waves + buoyancy heave. Pass: crest-trough ≈ 2A; equilibrium draft ≈ m/(ρA).
- **E4 Pipe network / fan laws.** Open `pipenetwork` (and `fanlaws`) → set flows/curves. Expect: node pressures / scaled fan curve. Pass: continuity holds; fan power ∝ speed³.
- **E5 Gas dynamics.** Open `gasdynamics` → set Mach, area ratio. Expect: isentropic ratios / normal-shock jump. Pass: matches isentropic tables.

## F. Thermal / HVAC
- **F1 HVAC load.** Open `hvac` → set room + ΔT. Expect: heating/cooling load. Pass: Q ∝ UAΔT.
- **F2 Heat pump.** Open `heatpump` → set source/sink temps. Expect: COP. Pass: COP ≤ Carnot, > 1.
- **F3 Thermal expansion.** Open `thermalexpansion` → set L, α, ΔT. Expect: ΔL. Pass: ΔL = αLΔT.

## G. Electrical / EM / signals
- **G1 Radar.** Open `sensors` → set Pt, gain, RCS, range. Expect: received power + SNR + detection. Pass: Pr ∝ R⁻⁴; detection toggles at the range limit.
- **G2 Transmission line.** Open `transmissionline` → set Z0, load. Expect: VSWR, reflection. Pass: matched load → VSWR 1.
- **G3 Filter + FFT.** Open `filter` then `fft` → set cutoff / signal. Expect: frequency response + spectrum. Pass: −3 dB at cutoff; FFT peak at the input tone.
- **G4 Three-phase / rectifier.** Open `threephase` (and `rectifier`) → set voltages. Expect: line/phase values; rectified ripple. Pass: V_line ≈ √3·V_phase.

## H. Bio & chem (the visualization-heavy section)
- **H1 Molecular viewer — representation + coloring.** Open genetics/biostruct → load a demo structure →
  set representation = `ribbon`, color = `Secondary structure`. Expect: a cartoon/ribbon colored helix-red/sheet-yellow/coil-grey. Pass: renders, colors correct.
- **H2 Trajectory playback.** Same viewer → attach demo trajectory → play. Expect: the structure animates across frames. Pass: frame slider advances, geometry moves.
- **H3 Molecular surface (SES).** Viewer → representation = surface, mode = `SES`, set probe radius. Expect: a smooth solvent-excluded surface. Pass: enclosing the vdW solid, re-entrant.
- **H4 Genetics tool.** Open `genetics` → pick a bio panel → run a small sequence task. Expect: a result (alignment/structure/etc.). Pass: result renders, no crash.
- **H5 PPI interactome.** Open `ppi` → set hosts/pathogens, analysis = degree centrality. Expect: a network graph + top-centrality readout. Pass: hub has max centrality.
- **H6 Reaction dynamics.** Open `reactdyn` → set preset + method → run. Expect: a reacting MD trajectory + energy plot. Pass: energy conserved (NVE), no blow-up.
- **H7 Enzyme kinetics.** Open `enzymekinetics` → set Vmax, Km, [S]. Expect: Michaelis–Menten curve. Pass: v → Vmax as [S]→∞.

## I. Sim & analysis (the cross-cutting enablers)
- **I1 UQ.** Open `uq` → set 2 input distributions, model, samples=1000 → run. Expect: output histogram + Sobol bars + FORM β/Pf. Pass: larger-variance input → larger Sobol index.
- **I2 ROM.** Open `rom` → set snapshots/rank → run. Expect: POD modes + error-vs-rank. Pass: error falls with rank.
- **I3 Co-simulation.** Open `cosim` → set H, scheme=Gauss-Seidel → run. Expect: coupled signal time-series + residual. Pass: residual → 0; tracks the monolithic ref at small H.
- **I4 Mission-sim.** Open `missionsim` → set blue/red, Pk, Lanchester. Expect: plan view + force-vs-time. Pass: Lanchester conserves a·A²−b·B².
- **I5 Autonomy V&V.** Open `autonomy` → set a scenario → run. Expect: pass/fail margins + clearance. Pass: margin sign matches the geometry.
- **I6 Photogrammetry.** Open `photogrammetry` → set cameras/points/noise → run. Expect: recovered cloud + reproj error. Pass: noise-free recovers to ~0 px.

## J. Machine design
- **J1 Gear train.** Open `gears` → set teeth/module. Expect: ratio + geometry. Pass: ratio = N2/N1.
- **J2 Spring.** Open `springs` → set wire/coil. Expect: rate k + a 3-D spring solid. Pass: k = Gd⁴/(8D³n); renders.
- **J3 Cam dynamics.** Open `camdynamics` → set profile/speed. Expect: lift/velocity/accel curves. Pass: continuous, no spikes.

## K. Neuro
- **K1 Stimulation.** Open `neuro` → set electrode/pulse → strength-duration. Expect: chronaxie/rheobase curve. Pass: matches the SD relation.

## L. Cross-cutting — the AI-drivability itself
- **L1 Open everything.** `OpenWorkbench{id}` for every registered id in turn. Pass: every panel opens, none crash.
- **L2 Discover + drive.** `ListCommands` + `ListControls{workbench}` on 5 panels → set + run by the listed names. Pass: every advertised name is settable; run succeeds.
- **L3 Viewport.** `set_view{iso}` → `add_sketch_point` ×4 → `extrude_sketch`. Pass: camera snaps; a part appears — all without a mouse.
- **L4 Bad input is safe.** `SetControl` with a junk value / unknown name on several panels. Pass: in-panel warning, never a panic.

---

### Standing gaps to expect (early adds)
1. **Result read-back** — a bridge command to return a panel's computed readout to the feed (so pass/fail auto-checks).
2. **3-D entity select-by-id** — no pick model yet (flagged open in the viewport work).
3. **Structure/file loading by path** — bridge currently uses synthetic demos; loading a real PDB/STL by path needs a path-arg command (file dialogs don't drive headless).
4. **Deeply-nested sub-forms** — meshtoolbox sub-panels, genetics' 15 sub-tools, astro per-stage — settable surface still to wire.
