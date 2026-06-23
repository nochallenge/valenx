# valenx Capability Roadmap — New Tracks (2026-06-22)

Four new capability tracks, sequenced by impact ÷ effort, all built to one principle.

## Guiding principle: all in-house, port-permissive

valenx is `MIT OR Apache-2.0`. Default to **in-house Rust**. When a reference implementation
exists, the action depends *only* on its license:

| Source license | Action | Examples |
|---|---|---|
| **Permissive** (BSD / MIT / Apache / MPL-2.0) | **Port the code directly into Rust**, keep the original copyright + license text in `THIRD-PARTY-NOTICES` (attribution is the only obligation). | COLMAP (BSD-3), NREL CCBlade (Apache-2.0), BYU FLOW-Lab CCBlade.jl/VortexLattice.jl (MIT), `peng_quad` (MIT) |
| **Copyleft** (GPL / LGPL / non-commercial) | **Do NOT port the code.** Reimplement from the published *math* (clean-room) — still in-house — or use as an optional, clearly-labeled arm's-length subprocess adapter. | Blender (GPL), XFOIL/AVL/QPROP/QBlade (GPL/NC), JSBSim (LGPL), FAST-UAV (GPL) |

This keeps valenx permissively licensed and genuinely in-house while still standing on the best
prior art. A line-by-line Rust port of permissive code is fine; a port of GPL code is a derivative
work and is not.

---

## Track A — Black-Hole / Relativity workbench  ·  **fastest payoff, do first**

**Status:** the engine is *done and validated* — `valenx-relativity` (#868) ships Schwarzschild /
Kerr / Reissner-Nordström / Kerr-Newman, autodiff curvature tensors, geodesic integration,
horizons / ergosphere / ISCO / photon-sphere / redshift observables, a **shadow + lensing renderer**
(`render_shadow → ShadowImage`), and Hawking thermodynamics. There is **no dedicated workbench**
surfacing it. This track is almost entirely UI + viewport/image wiring — low risk, high "wow".

- **P1 — Observables panel:** spacetime picker + (M, spin a, charge Q); display horizons, ergosphere,
  ISCO, photon sphere, shadow radius, redshift, Hawking T / entropy / luminosity / lifetime.
- **P2 — Shadow / lensing image:** `render_shadow → ShadowImage` into the image viewer (reuse the
  path-tracer render panel); Sgr A* / M87* presets.
- **P3 — Geodesics in the 3-D viewport:** trace photon + massive orbits (`integrate_geodesic`); show
  precession, light deflection.
- **P4 — extras:** accretion-disk Doppler/redshift coloring; a "clock at r vs ∞" time-dilation tool.
- **AI-drivable:** `BlackHole { spacetime, mass, spin, charge }` agent command + named widgets.

**Effort:** low. **Impact:** high (showcase). **Dependencies:** none (engine exists).

---

## Track B — Meshgen (text → 3D)  ·  **already planned**

**Status:** design spec + Phase-1 plan committed
(`docs/design/2026-06-22-meshgen-byo-llm-design.md` + `.../plans/2026-06-22-meshgen-byo-llm-plan.md`).
BYO mesh-LLM; in-house runtime behind a `MeshLlm` trait (llama-cpp-2 first → candle pure-Rust),
quantized-OBJ ↔ `valenx_mesh::Mesh` codec, "Text → 3D" workbench feeding the editor. Feature-gated,
off by default. **Ready to build.**

---

## Track C — Drone: design → simulate  ·  **in-house, port permissive aero**

**Status:** feasibility brief complete. valenx already owns **~70% of the spine in-house**:
`valenx-astro` (6-DOF rigid-body + quaternion attitude + PD control in `flight6dof.rs`, Mach drag,
propulsion), `valenx-vehicle` (battery/energy), `valenx-aero` (immersed-boundary Navier-Stokes
CFD), `valenx-mbd` (planar multibody). The **one gap is aerodynamics** — no BEMT / vortex-lattice /
airfoil code anywhere. That is the headline new capability, and the best references are permissive,
so it's a clean port, not clean-room.

- **P1 — `valenx-rotor` (BEMT):** blade-element-momentum propeller/rotor model, ported from NREL
  CCBlade (Apache) / BYU CCBlade.jl (MIT). Prandtl tip/hub loss + Glauert/Buhl high-thrust + a
  guaranteed-convergence 1-D residual. Thrust/torque/power/efficiency. *No Rust BEMT exists — this is
  the novel piece.*
- **P2 — multirotor force aggregator** on the existing `flight6dof` core (per-rotor thrust/torque ×
  moment arms) → full 6-DOF multirotor flight sim.
- **P3 — fixed-wing:** clean-room **VLM** (port BYU VortexLattice.jl, MIT) → lift, induced drag,
  stability derivatives; **tabulated airfoil polars** (XFOIL only as an optional GPL subprocess —
  never linked).
- **P4 — endurance + parametric airframe:** momentum-theory hover power + Peukert on `valenx-vehicle`;
  airframe via the existing CAD/curves/assembly crates.
- **P5 — generative auto-design:** NSGA-II / graph-grammar topology+controller search on
  `valenx-optimize`, scored by the sim. Methods from MIT/TU Delft papers (no permissive code exists —
  clean-room, still in-house).
- **Validation:** cross-check vs `peng_quad` (MIT) + `valenx-aero` CFD. **Optional adapter:** PX4 /
  Gazebo SITL over the permissive `mavlink` crate (arm's-length).

**Effort:** medium-high (`valenx-rotor` is real work). **Impact:** high (a whole new domain).

---

## Track D — Phone → 3-D scan  ·  **in-house photogrammetry, port COLMAP (BSD)**

**Status:** architecture set; crate-level detail finalizes with the scanning research brief. Three
parts, only the middle is hard:

1. **Capture + transport (mostly built):** extend the existing `valenx-remote` LAN+PIN web bridge —
   the phone browser captures camera via `getUserMedia` (frames + `MediaRecorder` video) and POSTs to
   the desktop; plus a file-upload for a pre-recorded video scan. (Same pipe, reversed direction.)
2. **Reconstruction (the work) — in-house Rust, porting COLMAP (BSD-3):**
   - **P1 — SfM front-end:** features (ORB/SIFT) → matching → essential-matrix pose → triangulation →
     **sparse point cloud**. Use Rust numerics (nalgebra/faer) for bundle-adjustment least-squares
     instead of vendoring Ceres.
   - **P2 — dense MVS + meshing:** dense stereo → point cloud → Poisson surface reconstruction (heavier).
   - **Optional 2nd backend:** **Gaussian splatting** via the Rust + wgpu `brush` engine (the
     "new-tech" path) behind the same interface — also in-house.
3. **Integration:** point cloud / mesh → the viewport + editor (point-cloud display likely needs adding;
   mesh import exists).

**Honest limits:** **web = camera frames/video only — no LiDAR** (that needs a native iOS app);
reconstruction is **GPU-heavy** and the hard part; rooms harder than objects; realistic UX is
"capture → process → result in seconds-to-minutes," not live "watch it build."

**Effort:** high (a real CV pipeline). **Impact:** high + distinctive.

---

## Sequencing recommendation

1. **Now, in parallel (low-risk):** **Track A black-hole** (fastest real payoff — engine's done) and
   **Track B meshgen** (already planned).
2. **Next:** **Track C drones** — clear brief; `valenx-rotor` (BEMT) is the headline build.
3. **Then:** **Track D scanning** — biggest CV undertaking; phase the COLMAP port (SfM first).

**Cross-cutting (every track):**
- **AI-drivable-first** — each workbench ships its agent-bridge command + named widgets.
- The **workbench-UX queue** (close/collapse/pop-out/dock/tabs, held local) continues underneath.
- Each new crate: feature-gate heavy/optional deps; `THIRD-PARTY-NOTICES` for any ported permissive code.

## Per-feature plan status
- **Meshgen** — plan written: `docs/design/plans/2026-06-22-meshgen-byo-llm-plan.md` ✓
- **Black-hole** — spec + plan: *next to write* (fastest payoff).
- **Drones** — spec + plan: to write (brief in hand).
- **Scanning** — spec + plan: to write once the research brief lands.
