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

## Track D — Phone → 3-D scan  ·  **in-house Rust port of COLMAP's SfM (+ brush splatting; dense-MVS the frontier)**

**Directive (user, 2026-06-22):** rewrite COLMAP in Rust. It's **BSD-3-Clause (permissive)**, so the
port is legitimate — keep COLMAP's copyright + BSD text in `THIRD-PARTY-NOTICES`.

**Honest scope:** COLMAP is a *pipeline*, not one program. The **sparse SfM front-end** ports cleanly
and is the high-value core that powers the scan. The **dense MVS stage** (PatchMatch stereo + fusion)
is the one piece with no Rust foundation — a genuine months-long from-scratch effort, and where COLMAP
relies on CUDA. Plan: **port the SfM core now**; get photoreal output from the in-house Rust **`brush`**
splatter (fed by our own poses) instead of porting CUDA dense-MVS; leave true dense-MVS meshing as a
later frontier (or an optional COLMAP adapter for users who must have it).

**Gating item:** the capture page must be served over **HTTPS** — `getUserMedia` needs a secure context,
so the current plain-HTTP `valenx-remote` bridge won't get camera access (add a self-signed cert trusted
on the phone). LiDAR stays web-impossible (native iOS only).

**Capture + transport:** extend `valenx-remote` over HTTPS — `getUserMedia` + `<video>→<canvas>→toBlob`
frames + `MediaRecorder` video (MP4 iOS / WebM Android) POSTed to the desktop, plus a file-upload for a
pre-recorded scan. Reconstruction is an **async job** ("capture → upload → reconstruct → notify").

**Build — new crate `valenx-photogrammetry`, COLMAP's SfM ported stage by stage:**
- **D0 — Instant depth preview (in-house, days):** `candle` + Depth-Anything-V2 Small (Apache) → a rough
  2.5-D point cloud while the full SfM runs. Quick win; proves the phone→desktop loop.
- **D1 — feature extraction** (SIFT/ORB — port COLMAP's stage or build on `kornia-rs`).
- **D2 — matching** (descriptor matching + RANSAC geometric verification).
- **D3 — two-view geometry** (essential/fundamental matrix → relative pose → initial triangulation).
- **D4 — incremental mapper** (PnP registration + triangulation + incremental loop) → camera poses +
  a sparse point cloud.
- **D5 — bundle adjustment** (Levenberg-Marquardt on nalgebra/faer sparse, or `apex-solver`) — the
  riskiest link; jointly refines poses + points.
- **D6 — photoreal output:** feed the SfM poses to the in-house Rust **`brush`** Gaussian-splat trainer
  (Apache, wgpu) → a photoreal model. (`brush` pins forked wgpu/egui — `cargo deny` + reconcile.)
- **Frontier (optional, later):** dense MVS (port COLMAP's PatchMatch stereo) for true dense meshes —
  the months-long piece — or an optional external COLMAP adapter.

**Integration:** sparse cloud / splat / mesh → the viewport + editor (point-cloud + splat display need
adding; mesh import exists).

**Avoid (non-commercial / copyleft):** Inria 3DGS + SuGaR (non-commercial), DA-V2 Base/Large (CC-BY-NC),
OpenMVS (AGPL — optional user-installed adapter only, never bundled).

**Honest scale + limits:** the **largest of the four tracks** — weeks for the SfM core, dense-MVS its own
project; objects scan far better than rooms; reconstruction is async; no LiDAR over the web.

---

## Sequencing recommendation

1. **Now, in parallel (low-risk):** **Track A black-hole** (fastest real payoff — engine's done) and
   **Track B meshgen** (already planned).
2. **Next:** **Track C drones** — clear brief; `valenx-rotor` (BEMT) is the headline build.
3. **Then:** **Track D scanning** — in-house Rust port of COLMAP's SfM (new `valenx-photogrammetry`, phased) + `brush` photoreal; the largest track. HTTPS for capture first.

**Cross-cutting (every track):**
- **AI-drivable-first** — each workbench ships its agent-bridge command + named widgets.
- The **workbench-UX queue** (close/collapse/pop-out/dock/tabs, held local) continues underneath.
- Each new crate: feature-gate heavy/optional deps; `THIRD-PARTY-NOTICES` for any ported permissive code.

## Per-feature plan status
- **Meshgen** — plan written: `docs/design/plans/2026-06-22-meshgen-byo-llm-plan.md` ✓
- **Black-hole** — plan written: `docs/design/plans/2026-06-22-blackhole-workbench-plan.md` ✓
- **Drones** — research brief in hand (Track C is the phased plan); full TDD plan to write on demand.
- **Scanning** — in-house COLMAP-SfM port (Track D, new `valenx-photogrammetry` crate); full TDD plan to write next.
