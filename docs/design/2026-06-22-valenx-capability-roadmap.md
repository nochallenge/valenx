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

## Track D — Phone → 3-D scan  ·  **tiered: in-house depth + splatting, COLMAP adapter for meshes**

**Status:** research complete (full brief in the session log). **Correction to the first framing:** a
*full* in-house Rust photogrammetry-to-dense-mesh pipeline is **not feasible today** — the dense
multi-view-stereo stage has **zero Rust implementations**, so a wholesale COLMAP port would be a
multi-year from-scratch CV effort. The genuinely in-house path is **monocular depth + Gaussian
splatting** (both Apache-licensed, pure Rust); **COLMAP stays an optional adapter** for users who want
a precise geometric mesh (BSD — portable in principle, impractical to port in full).

**Gating item:** the capture page must be served over **HTTPS** — `getUserMedia` requires a secure
context, so the current plain-HTTP `valenx-remote` bridge won't get camera access (add a self-signed
cert trusted on the phone). LiDAR stays web-impossible (native iOS only).

**Capture + transport:** extend `valenx-remote` over HTTPS — `getUserMedia` + `<video>→<canvas>→toBlob`
frames + `MediaRecorder` video (MP4 on iOS / WebM on Android) POSTed to the desktop, plus a file-upload
for a pre-recorded scan. Reconstruction runs as an **async job** ("capture → upload → reconstruct →
notify"), never a live viewfinder.

**Tiered build (ship in this order):**
- **D1 — Instant depth preview (in-house, fastest):** phone photo → `candle` + **Depth-Anything-V2
  Small** (Apache-2.0, BYO weights) → back-projected point cloud in the viewport. Honestly a 2.5-D,
  scale-ambiguous preview — not metrology. Proves the phone→desktop loop end-to-end. (~days)
- **D2 — Gaussian splatting (in-house flagship, the "wow"):** integrate **`brush`** (Apache-2.0,
  Rust + **wgpu** — valenx's exact GPU stack, CUDA-free, runs on AMD/Intel/Apple) to train a photoreal
  splat from a photo set / video. Needs posed input (COLMAP-sparse adapter, CPU, supplies poses) and
  outputs render-only `.ply` splats (no mesh yet). Caveat: `brush` pins git forks of wgpu/egui/naga —
  run `cargo deny` + reconcile against valenx's wgpu. (~weeks)
- **D3 — COLMAP adapter → true mesh (optional external):** arm's-length subprocess to user-installed
  **COLMAP (BSD)** for full SfM→dense→mesh; honest that **dense needs NVIDIA CUDA** (CPU-only users get
  sparse/poses only). Also the pose-provider for D2.

**Integration:** point cloud / splat / mesh → the viewport + editor (point-cloud + splat display need
adding; mesh import exists).

**Avoid (non-commercial / copyleft):** Inria 3DGS + SuGaR (non-commercial), DA-V2 Base/Large
(CC-BY-NC), OpenMVS (AGPL — only ever an optional user-installed adapter, never bundled).

**Honest limits:** async (seconds for depth; minutes–tens-of-minutes for splatting/COLMAP); **objects
scan far better than rooms** (rooms really want LiDAR, which the web can't give); no LiDAR over the web.

**Effort:** D1 low · D2 medium · D3 low (adapter). **Impact:** high + distinctive.

---

## Sequencing recommendation

1. **Now, in parallel (low-risk):** **Track A black-hole** (fastest real payoff — engine's done) and
   **Track B meshgen** (already planned).
2. **Next:** **Track C drones** — clear brief; `valenx-rotor` (BEMT) is the headline build.
3. **Then:** **Track D scanning** — tiered: depth-preview (in-house) → `brush` splatting (in-house) → COLMAP adapter (real meshes). Serve the capture page over HTTPS first.

**Cross-cutting (every track):**
- **AI-drivable-first** — each workbench ships its agent-bridge command + named widgets.
- The **workbench-UX queue** (close/collapse/pop-out/dock/tabs, held local) continues underneath.
- Each new crate: feature-gate heavy/optional deps; `THIRD-PARTY-NOTICES` for any ported permissive code.

## Per-feature plan status
- **Meshgen** — plan written: `docs/design/plans/2026-06-22-meshgen-byo-llm-plan.md` ✓
- **Black-hole** — plan written: `docs/design/plans/2026-06-22-blackhole-workbench-plan.md` ✓
- **Drones** — research brief in hand (Track C is the phased plan); full TDD plan to write on demand.
- **Scanning** — research brief complete; tiered plan in Track D. D1 (depth-preview) TDD plan to write next.
