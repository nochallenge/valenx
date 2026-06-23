# valenx-meshgen — Bring-Your-Own Mesh-LLM (Text ↔ 3D) — Design

**Status:** Approved design (2026-06-22).
**Implementation plan:** `docs/design/plans/2026-06-22-meshgen-byo-llm-plan.md` (Phase 1).

**Goal:** Let a user generate *editable* 3D meshes from a text prompt — and describe / refine
existing meshes — using a local mesh-generation LLM **they** download and plug in. valenx ships
the entire pipeline; it never ships the weights.

**Architecture (one line):** A feature-gated `valenx-meshgen` crate exposes a `MeshLlm` trait with
swappable backends (llama-cpp-2 first, candle later); a codec turns the model's quantized-OBJ token
stream into a *validated* valenx `Mesh`; a "Text → 3D" workbench runs generation on a background
thread and drops the result into the existing viewport/editor.

**Tech stack:** Rust · llama-cpp-2 (GGUF inference) → candle (pure Rust) behind a trait · nalgebra ·
`valenx-mesh` (Mesh + validate/weld/repair) · eframe/egui (workbench) · the existing background-job
and agent-bridge infrastructure.

---

## 1. Motivation

valenx is a precision engineering/CAD tool with a real geometry kernel and a 3-D editor. What it
lacks is a *fast on-ramp to a starting shape*. Mesh-generation LLMs (LLaMA-Mesh and successors)
turn a text prompt into clean, low-poly, **editable topology** — the opposite of the messy soup that
diffusion/NeRF text-to-3D produces. That is exactly the right input for a precision editor: the model
gives you a rough draft, and valenx's existing push/pull/gizmo/CSG tools make it precise.

We deliberately do **not** bundle or run a model we ship. An 8B model is ~5 GB and GPU-hungry, and the
weights carry their own licenses. Instead valenx ships the *plumbing* and the user supplies the model
file (BYO). This keeps valenx small, license-clean, and forward-compatible with better models as they
appear.

## 2. Non-goals (YAGNI)

- **No bundled weights, no training/fine-tuning.** valenx is the runtime + UI only.
- **Not a general chat assistant.** Scoped to mesh generation / description / refinement.
- **Not high-poly or photoreal.** Output is low-poly concept topology *by design* — it feeds the editor.
- **Not a cloud service.** Inference is local, on the user's machine, with their model.
- `describe` and `chat-refine` are later phases (P3), not v1.

## 3. Architecture & components

### 3.1 `valenx-meshgen` crate (new, feature-gated)

A single crate owning the runtime abstraction, the backends, the codec, and the model profile. The
whole crate — and the heavy `llama-cpp-2` dependency — lives behind a cargo feature `meshgen`,
**off by default**, so the normal workspace build and CI never compile the LLM stack.

### 3.2 The `MeshLlm` trait — the runtime seam

```rust
/// A loadable mesh-generation language model. Backends are swappable; the app
/// never knows which one it holds.
pub trait MeshLlm: Send {
    /// Load a model from disk onto `device`.
    fn load(config: &ModelConfig, device: Device) -> Result<Self, MeshGenError>
    where
        Self: Sized;

    /// Stream a completion for `prompt`. `on_token` is called for each decoded
    /// token; returning `ControlFlow::Break` cancels generation promptly.
    fn generate(
        &mut self,
        prompt: &str,
        params: &GenParams,
        on_token: &mut dyn FnMut(&str) -> std::ops::ControlFlow<()>,
    ) -> Result<(), MeshGenError>;
}
```

- **`LlamaCppBackend`** (Phase 1) — wraps `llama-cpp-2`; loads a GGUF; streams tokens. Runs on CPU or
  any GPU llama.cpp supports (CUDA / Vulkan / Metal).
- **`CandleBackend`** (Phase 4) — pure-Rust (HuggingFace `candle`); safetensors; the "in-house" target.
- **`MockBackend`** (test-only) — replays a canned token stream, so the codec and the app state machine
  are testable with no model present.

### 3.3 The codec — quantized OBJ ⇄ valenx `Mesh`

A `codec` module with two directions:

- **`decode`** — parse the model's streamed `v x y z` / `f a b c` text, dequantize the integer grid
  (per the model profile) into a valenx `Mesh`. This is essentially the existing OBJ parse plus a
  dequant step. It is **streaming-tolerant**: it accepts partial input and *skips* malformed/incomplete
  lines rather than panicking (consistent with the parser hardening already in the codebase).
- **`encode`** — serialize a `Mesh` to quantized-OBJ text, for the `describe` direction (P3).

### 3.4 `ModelProfile` / `ModelConfig`

- **`ModelProfile`** — prompt template, quantization resolution (e.g. a 0–127 coordinate grid),
  chat/turn format, special tokens. LLaMA-Mesh's defaults are baked in as a constant; the struct is the
  seam any other mesh-LLM plugs into.
- **`ModelConfig`** — path to the weights file + the `ModelProfile` + device selection.

### 3.5 App integration (`valenx-app`)

- A `meshgen_workbench.rs` panel: prompt box, **Generate**, a cancellable streaming progress indicator,
  and settings (model path, backend, device).
- Generation runs on a **background thread** using the existing run/job infrastructure; tokens stream
  back to the UI. On completion the decoded `Mesh` is run through `valenx-mesh`'s **validate / weld /
  repair** (LLM mesh output is often non-manifold) and loaded into the viewport via the existing
  `mesh_loader` path — where the editor takes over.
- One or more **agent-bridge commands** (`MeshGen { prompt }`) so the feature is AI-drivable, per the
  standing AI-drivable-first priority.

## 4. Data flow

**Generate:** `prompt → workbench → background thread → MeshLlm.generate → token stream → codec.decode
(incremental) → Mesh → validate/weld/repair → viewport → editor`.

**Describe (P3):** `viewport Mesh → codec.encode → prompt → MeshLlm.generate → text → panel`.

## 5. Error handling & honesty

- **No model configured →** the workbench shows a clear "download a mesh LLM and set its path" state.
  Never a placeholder mesh. (BLOCKED-not-faked contract.)
- **Load failure** (bad path, wrong format, OOM) → a typed `MeshGenError` surfaced in the panel.
- **Invalid generated geometry →** run through validate/weld/repair; if still degenerate, report it
  honestly ("the model produced an invalid mesh — try again or adjust the prompt"). Never silently show
  garbage.
- **Cancellation** → `on_token` returns `ControlFlow::Break`; generation stops promptly and the thread
  unwinds cleanly.

## 6. Feature-gating & build

- `valenx-meshgen` + `llama-cpp-2` are behind the `meshgen` cargo feature, **off by default**. The
  default workspace build and the existing CI matrix never compile the LLM stack.
- A separate, opt-in CI job builds `--features meshgen` (CPU backend) so the code stays green without
  burdening every run.
- **Build risk:** `llama-cpp-2` compiles llama.cpp (CUDA/Vulkan toolchain friction on Windows). The
  feature gate contains it, and Phase 1 proves the codec + a CPU build **headless** before any GUI work.

## 7. Testing strategy

- **Codec:** encode→decode round-trip; golden quantized-OBJ fixtures; malformed/partial-stream input
  (asserts no panic, graceful skip); dequant accuracy bounds.
- **MeshLlm:** `MockBackend` replays a canned token stream → fully deterministic tests with no weights.
- **Headless integration (P1):** with a tiny real GGUF supplied locally (BYO; **not** in CI),
  `generate("a cube") → Mesh` smoke test, `#[ignore]`d by default and run via an env-var gate.
- **Workbench (P2):** a headless-UI-style test of the panel state machine
  (no-model → loaded → generating → done → error).

## 8. Blender parallel track (clean-room Rust editor ops)

Separate from meshgen, noted because the generated meshes feed the editor. The right way to get
"Blender-grade" editing is to **extend `valenx-blender-mesh-ops`** — already a clean-room Rust
reimplementation of Blender-*style* ops (extrude region, bevel, inset, loop cut, bridge, boolean,
solidify) — with the high-value missing operations (a fuller bevel, remesh, etc.).

**Legal constraint (firm):** these are implemented from the *algorithms* (papers / documented behavior),
**never ported from Blender's GPL source.** A line-by-line translation of GPL C++ into Rust is a
derivative work and would force valenx to GPL — incompatible with its MIT/Apache license. Idea =
reimplementable; expression = not. This track gets its own spec/plan when prioritized; it is listed here
only so the roadmap is coherent (generate → edit with strong, license-clean ops).

## 9. Phasing

- **P1 — `valenx-meshgen` crate:** `MeshLlm` trait + `LlamaCppBackend` + codec + `MockBackend` +
  headless prompt→Mesh test. *(Covered by this plan.)*
- **P2 — Text→3D workbench:** background generation, streaming UI, validation pipeline, viewport load,
  settings, agent-bridge command.
- **P3 — describe + chat-refine:** mesh→text encode; conversational iterative edits.
- **P4 — candle backend:** pure-Rust backend behind the trait.

## 10. Risks

- **llama-cpp-2 build friction** — mitigated by the feature gate, a CPU-first path, and a headless P1.
- **Variable model quality** (BYO) — mitigated by the honest UX and the validate/weld/repair pipeline.
- **Streaming decode of partial meshes** — the codec must tolerate incomplete face lists mid-stream
  (covered by the partial-stream tests).
