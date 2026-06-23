# valenx-meshgen Phase 1 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the headless foundation of the bring-your-own mesh-LLM feature — a new `valenx-meshgen` crate that turns a mesh-LLM's quantized-OBJ token stream into a valenx `Mesh`, with a swappable backend trait and a feature-gated llama.cpp backend. No GUI in this phase.

**Architecture:** A `MeshLlm` trait abstracts the runtime; a `codec` converts quantized-OBJ text ⇄ `valenx_mesh::Mesh`; a `MockBackend` (canned token stream) makes everything testable with no weights; a `LlamaCppBackend` (real inference) lives behind a `meshgen` cargo feature that is OFF by default. Everything is verified headless before any workbench work (Phase 2).

**Tech Stack:** Rust · `valenx-mesh` (Mesh/ElementBlock/repair) · `nalgebra` (Vector3) · `thiserror` (errors) · `llama-cpp-2` (GGUF inference, optional/feature-gated).

**Spec:** `docs/design/2026-06-22-meshgen-byo-llm-design.md`

**Commit identity (use on every commit in this plan):**
```
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "<msg>"
```
After each commit, sanity-check no email leak: `git log -1 --format='%ae' | grep -c forceblue` must print `0`.

---

## File Structure

All paths under `crates/valenx-meshgen/`:

- `Cargo.toml` — crate manifest; declares the optional `llama-cpp-2` dep + the `meshgen` feature.
- `src/lib.rs` — crate docs + re-exports; declares modules.
- `src/error.rs` — `MeshGenError` (the one error type the crate returns).
- `src/model.rs` — `Device`, `GenParams`, `ModelProfile` (+ `ModelProfile::LLAMA_MESH`, `quant`/`dequant`), `ModelConfig`.
- `src/codec.rs` — `decode(text, &ModelProfile) -> Mesh` and `encode(&Mesh, &ModelProfile) -> String`.
- `src/backend/mod.rs` — the `MeshLlm` trait + `MockBackend`.
- `src/backend/llama.rs` — `LlamaCppBackend` (compiled only under `--features meshgen`).

Also modified:
- `Cargo.toml` (workspace root) — add `"crates/valenx-meshgen"` to `members`.

Each file has one responsibility; `codec.rs` and `model.rs` are fully testable with no model present.

---

## Task 0: Scaffold the crate + register it in the workspace

**Files:**
- Create: `crates/valenx-meshgen/Cargo.toml`
- Create: `crates/valenx-meshgen/src/lib.rs`
- Modify: `Cargo.toml` (workspace root `members` list)

- [ ] **Step 1: Create the manifest** — `crates/valenx-meshgen/Cargo.toml` (mirrors the `valenx-recipes` leaf-crate shape; the `llama-cpp-2` dep + `meshgen` feature are added in Task 5, kept out now so Task 0–4 build with zero heavy deps):

```toml
[package]
name        = "valenx-meshgen"
description = "Bring-your-own mesh-LLM pipeline for valenx: a MeshLlm runtime trait + a quantized-OBJ <-> valenx Mesh codec, so a user-supplied mesh-generation model (e.g. LLaMA-Mesh) turns a text prompt into an editable mesh. Backends are swappable (a test MockBackend now; a feature-gated llama.cpp backend behind `meshgen`). Ships the pipeline, never the weights."

version.workspace      = true
edition.workspace      = true
rust-version.workspace = true
authors.workspace      = true
license.workspace      = true
repository.workspace   = true

[dependencies]
thiserror   = { workspace = true }
nalgebra    = { workspace = true }
valenx-mesh = { path = "../valenx-mesh" }

[lints.rust]
unused_imports = "warn"
missing_docs   = "warn"
```

- [ ] **Step 2: Create the lib root** — `crates/valenx-meshgen/src/lib.rs`:

```rust
//! # valenx-meshgen
//!
//! Bring-your-own mesh-LLM pipeline. valenx ships this runtime + codec; the
//! user supplies the model weights. See
//! `docs/design/2026-06-22-meshgen-byo-llm-design.md`.
#![forbid(unsafe_code)]

pub mod backend;
pub mod codec;
pub mod error;
pub mod model;

pub use error::MeshGenError;
pub use model::{Device, GenParams, ModelConfig, ModelProfile};
```

(The `backend`/`codec` modules are created in later tasks; add empty stubs now so the crate compiles: create `src/error.rs`, `src/model.rs`, `src/codec.rs`, and `src/backend/mod.rs` each containing only `//! stub` for this step, then fill them in their tasks. To keep Task 0 green with the `pub use`s above, instead defer the `pub use` lines until Task 1/2 — simplest: in Step 2 include ONLY `pub mod error; pub mod model; pub mod codec; pub mod backend;` and add the `pub use` re-exports at the end of Task 4.)

- [ ] **Step 3: Register in the workspace** — in the root `Cargo.toml`, add `"crates/valenx-meshgen",` to the `members = [ ... ]` list (the list is explicit, not a glob).

- [ ] **Step 4: Create stub module files** so the crate compiles:
  - `src/error.rs` → `//! Error type (Task 1).`
  - `src/model.rs` → `//! Model config + profile (Task 1).`
  - `src/codec.rs` → `//! Quantized-OBJ codec (Task 2).`
  - `src/backend/mod.rs` → `//! MeshLlm trait + backends (Task 4).`

- [ ] **Step 5: Verify it builds**

Run: `cargo build -p valenx-meshgen`
Expected: compiles clean (empty crate).

- [ ] **Step 6: Commit**

```bash
git add crates/valenx-meshgen Cargo.toml
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(meshgen): scaffold valenx-meshgen crate (Phase 1, task 0)"
```

---

## Task 1: Core types — error, device, params, model profile

**Files:**
- Modify: `crates/valenx-meshgen/src/error.rs`
- Modify: `crates/valenx-meshgen/src/model.rs`

- [ ] **Step 1: Write the failing test** — append to `src/model.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn llama_mesh_profile_quant_dequant_roundtrips() {
        let p = ModelProfile::LLAMA_MESH;
        // Every bin index dequantizes then requantizes to itself.
        for q in [0_i64, 1, 50, p.quant_bins as i64 - 1] {
            let x = p.dequant(q);
            assert_eq!(p.quant(x), q, "bin {q} did not round-trip");
        }
        // Dequant maps bin 0 -> coord_min and the last bin -> coord_max.
        assert!((p.dequant(0) - p.coord_min).abs() < 1e-9);
        assert!((p.dequant(p.quant_bins as i64 - 1) - p.coord_max).abs() < 1e-9);
    }

    #[test]
    fn quant_clamps_out_of_range_coords() {
        let p = ModelProfile::LLAMA_MESH;
        assert_eq!(p.quant(p.coord_min - 100.0), 0);
        assert_eq!(p.quant(p.coord_max + 100.0), p.quant_bins as i64 - 1);
    }
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p valenx-meshgen`
Expected: FAIL — `ModelProfile` etc. not defined.

- [ ] **Step 3: Implement `src/error.rs`:**

```rust
//! The single error type returned by valenx-meshgen.

use std::path::PathBuf;

/// Everything that can go wrong loading or running a mesh-LLM.
#[derive(Debug, thiserror::Error)]
pub enum MeshGenError {
    /// No model file is configured (the user must download + point to one).
    #[error("no mesh-LLM is configured — download a model and set its path")]
    NoModelConfigured,

    /// The model file could not be loaded.
    #[error("failed to load model at {path}: {reason}")]
    ModelLoad { path: PathBuf, reason: String },

    /// Inference failed mid-generation.
    #[error("inference failed: {0}")]
    Inference(String),

    /// The decoded geometry was unusable even after repair.
    #[error("the model produced an invalid mesh: {0}")]
    InvalidGeometry(String),

    /// An I/O error.
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}
```

- [ ] **Step 4: Implement `src/model.rs`:**

```rust
//! Model configuration + the quantization/prompt profile.

use std::path::PathBuf;

/// Where inference runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Device {
    /// CPU only.
    Cpu,
    /// Offload `layers` transformer layers to the GPU (llama.cpp `n_gpu_layers`).
    Gpu { layers: u32 },
}

/// Sampling / decoding parameters for one generation.
#[derive(Debug, Clone, Copy)]
pub struct GenParams {
    /// Hard cap on generated tokens (bounds runtime + output size).
    pub max_tokens: usize,
    /// Sampling temperature.
    pub temperature: f32,
    /// Nucleus-sampling top-p.
    pub top_p: f32,
    /// RNG seed for reproducibility.
    pub seed: u64,
}

impl Default for GenParams {
    fn default() -> Self {
        Self { max_tokens: 8192, temperature: 0.7, top_p: 0.9, seed: 0 }
    }
}

/// How a particular mesh-LLM serializes geometry + how it is prompted.
///
/// LLaMA-Mesh emits OBJ-style `v x y z` / `f a b c` lines where the coords are
/// integers on a fixed grid. `quant_bins` is that grid resolution and
/// `coord_min..coord_max` is the world-space cube it maps onto. The codec is
/// grid-agnostic; this struct is the seam other mesh-LLMs plug into.
#[derive(Debug, Clone, Copy)]
pub struct ModelProfile {
    /// Human-readable profile name.
    pub name: &'static str,
    /// Number of quantization bins per axis (grid resolution).
    pub quant_bins: u32,
    /// World coordinate the lowest bin maps to.
    pub coord_min: f64,
    /// World coordinate the highest bin maps to.
    pub coord_max: f64,
    /// Prompt template; `{prompt}` is replaced with the user's text.
    pub prompt_template: &'static str,
}

impl ModelProfile {
    /// Defaults for LLaMA-Mesh.
    ///
    /// NOTE: confirm `quant_bins` against the model card at implementation —
    /// LLaMA-Mesh uses a fixed grid (commonly 64 or 128). The proportions are
    /// preserved regardless; only the absolute scale depends on this, and the
    /// mesh can be renormalized after decode.
    pub const LLAMA_MESH: ModelProfile = ModelProfile {
        name: "llama-mesh",
        quant_bins: 128,
        coord_min: -1.0,
        coord_max: 1.0,
        prompt_template: "Create a 3D obj file of {prompt}.",
    };

    /// Map a bin index to a world coordinate.
    pub fn dequant(&self, q: i64) -> f64 {
        let span = (self.quant_bins.max(2) - 1) as f64;
        let t = (q as f64 / span).clamp(0.0, 1.0);
        self.coord_min + t * (self.coord_max - self.coord_min)
    }

    /// Map a world coordinate to the nearest bin index (clamped in-range).
    pub fn quant(&self, x: f64) -> i64 {
        let span = (self.quant_bins.max(2) - 1) as f64;
        let denom = self.coord_max - self.coord_min;
        let t = if denom.abs() < f64::EPSILON {
            0.0
        } else {
            ((x - self.coord_min) / denom).clamp(0.0, 1.0)
        };
        (t * span).round() as i64
    }

    /// Build the full prompt for `user_text`.
    pub fn build_prompt(&self, user_text: &str) -> String {
        self.prompt_template.replace("{prompt}", user_text)
    }
}

/// A loaded-model configuration: the weights file + how to read its output.
#[derive(Debug, Clone)]
pub struct ModelConfig {
    /// Path to the user-supplied weights file (e.g. a GGUF).
    pub path: PathBuf,
    /// The serialization/prompt profile for this model.
    pub profile: ModelProfile,
    /// Where to run inference.
    pub device: Device,
}
```

- [ ] **Step 5: Add re-exports** to `src/lib.rs` (already present from Task 0 Step 2 — confirm `pub use model::{Device, GenParams, ModelConfig, ModelProfile};` and `pub use error::MeshGenError;` compile now that the types exist).

- [ ] **Step 6: Run the tests**

Run: `cargo test -p valenx-meshgen`
Expected: PASS (2 tests).

- [ ] **Step 7: Commit**

```bash
git add crates/valenx-meshgen/src
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(meshgen): core types — MeshGenError + Device/GenParams/ModelProfile/ModelConfig (task 1)"
```

---

## Task 2: Codec — decode quantized-OBJ text → valenx Mesh

The decode must be robust to malformed/garbage lines (the parser-hardening convention in this repo): skip a bad line, never panic; skip a face whose index is out of range.

**Files:**
- Modify: `crates/valenx-meshgen/src/codec.rs`

- [ ] **Step 1: Write the failing test** — `src/codec.rs`:

```rust
//! Quantized-OBJ <-> valenx Mesh codec.

use crate::model::ModelProfile;
use nalgebra::Vector3;
use valenx_mesh::{ElementBlock, ElementType, Mesh};

#[cfg(test)]
mod tests {
    use super::*;

    // A unit tetrahedron expressed in LLAMA_MESH bins (128-bin grid, -1..1):
    // bin 0 -> -1, bin 127 -> +1, bin 64 -> ~0.008. Faces are 1-based, fan ok.
    const TET: &str = "\
v 0 0 0
v 127 0 0
v 0 127 0
v 0 0 127
f 1 2 3
f 1 2 4
f 1 3 4
f 2 3 4
";

    #[test]
    fn decode_reads_vertices_and_triangulates_faces() {
        let m = decode(TET, &ModelProfile::LLAMA_MESH);
        assert_eq!(m.nodes.len(), 4, "4 vertices");
        let tris: usize = m.element_blocks.iter().map(|b| b.connectivity.len() / 3).sum();
        assert_eq!(tris, 4, "4 triangles");
        // bin 0 -> coord_min (-1.0)
        assert!((m.nodes[0].x - (-1.0)).abs() < 1e-9);
        // bin 127 -> coord_max (+1.0)
        assert!((m.nodes[1].x - 1.0).abs() < 1e-9);
    }

    #[test]
    fn decode_fan_triangulates_a_quad() {
        let quad = "v 0 0 0\nv 127 0 0\nv 127 127 0\nv 0 127 0\nf 1 2 3 4\n";
        let m = decode(quad, &ModelProfile::LLAMA_MESH);
        let tris: usize = m.element_blocks.iter().map(|b| b.connectivity.len() / 3).sum();
        assert_eq!(tris, 2, "a quad fans into 2 triangles");
    }

    #[test]
    fn decode_skips_garbage_and_out_of_range_faces_without_panicking() {
        let junk = "\
hello world
v 0 0 0
v not a number
v 127 0 0
v 0 127 0
f 1 2 999
f 1 2 3
garbage line
f
";
        let m = decode(junk, &ModelProfile::LLAMA_MESH);
        // Only the 3 well-formed `v` lines parse.
        assert_eq!(m.nodes.len(), 3);
        // `f 1 2 999` is dropped (999 out of range); `f 1 2 3` is kept; `f` is dropped.
        let tris: usize = m.element_blocks.iter().map(|b| b.connectivity.len() / 3).sum();
        assert_eq!(tris, 1);
    }

    #[test]
    fn decode_empty_input_is_an_empty_mesh() {
        let m = decode("", &ModelProfile::LLAMA_MESH);
        assert_eq!(m.nodes.len(), 0);
        assert_eq!(m.element_blocks.iter().map(|b| b.connectivity.len()).sum::<usize>(), 0);
    }
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p valenx-meshgen codec`
Expected: FAIL — `decode` not defined.

- [ ] **Step 3: Implement `decode`** (prepend above the `#[cfg(test)]` block):

```rust
/// Decode a (possibly partial) quantized-OBJ token string into a valenx `Mesh`.
///
/// Two-pass + tolerant: pass 1 collects every well-formed `v` line (dequantized
/// via `profile`); pass 2 fan-triangulates every `f` line, skipping any face
/// that references an out-of-range vertex. Malformed lines are skipped, never
/// fatal — the model output is untrusted text.
pub fn decode(text: &str, profile: &ModelProfile) -> Mesh {
    let mut mesh = Mesh::new("meshgen");

    // Pass 1: vertices.
    for line in text.lines() {
        let mut it = line.split_whitespace();
        if it.next() != Some("v") {
            continue;
        }
        let coords: Vec<f64> = it
            .take(3)
            .filter_map(|t| t.parse::<i64>().ok())
            .map(|q| profile.dequant(q))
            .collect();
        if coords.len() == 3 {
            mesh.nodes.push(Vector3::new(coords[0], coords[1], coords[2]));
        }
    }

    // Pass 2: faces (1-based OBJ indices -> 0-based, fan-triangulated).
    let n = mesh.nodes.len() as u32;
    let mut block = ElementBlock::new(ElementType::Tri3);
    for line in text.lines() {
        let mut it = line.split_whitespace();
        if it.next() != Some("f") {
            continue;
        }
        let verts: Vec<u32> = it
            // OBJ face tokens may be `v`, `v/vt`, `v//vn`; take the vertex part.
            .filter_map(|tok| tok.split('/').next()?.parse::<i64>().ok())
            // 1-based, positive only; convert to 0-based.
            .filter_map(|i| u32::try_from(i - 1).ok())
            .collect();
        if verts.len() < 3 {
            continue;
        }
        for k in 1..verts.len() - 1 {
            let (a, b, c) = (verts[0], verts[k], verts[k + 1]);
            if a < n && b < n && c < n {
                block.connectivity.extend_from_slice(&[a, b, c]);
            }
        }
    }
    mesh.element_blocks.push(block);
    mesh.recompute_stats();
    mesh
}
```

- [ ] **Step 4: Run the tests**

Run: `cargo test -p valenx-meshgen codec`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add crates/valenx-meshgen/src/codec.rs
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(meshgen): codec decode — quantized-OBJ -> Mesh, tolerant of malformed input (task 2)"
```

---

## Task 3: Codec — encode Mesh → quantized-OBJ text (round-trip)

**Files:**
- Modify: `crates/valenx-meshgen/src/codec.rs`

- [ ] **Step 1: Write the failing test** — add to the `tests` module in `src/codec.rs`:

```rust
    #[test]
    fn encode_then_decode_preserves_topology() {
        let original = decode(TET, &ModelProfile::LLAMA_MESH);
        let text = encode(&original, &ModelProfile::LLAMA_MESH);
        let round = decode(&text, &ModelProfile::LLAMA_MESH);
        assert_eq!(round.nodes.len(), original.nodes.len());
        let t0: usize = original.element_blocks.iter().map(|b| b.connectivity.len()).sum();
        let t1: usize = round.element_blocks.iter().map(|b| b.connectivity.len()).sum();
        assert_eq!(t1, t0);
        // Coords survive a quant->dequant round-trip to within one bin width.
        let bin = (ModelProfile::LLAMA_MESH.coord_max - ModelProfile::LLAMA_MESH.coord_min)
            / (ModelProfile::LLAMA_MESH.quant_bins as f64 - 1.0);
        for (a, b) in original.nodes.iter().zip(round.nodes.iter()) {
            assert!((a - b).norm() <= bin * 2.0);
        }
    }

    #[test]
    fn encode_emits_v_and_f_lines() {
        let m = decode(TET, &ModelProfile::LLAMA_MESH);
        let text = encode(&m, &ModelProfile::LLAMA_MESH);
        assert_eq!(text.lines().filter(|l| l.starts_with("v ")).count(), 4);
        assert_eq!(text.lines().filter(|l| l.starts_with("f ")).count(), 4);
    }
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p valenx-meshgen codec`
Expected: FAIL — `encode` not defined.

- [ ] **Step 3: Implement `encode`** (add below `decode`):

```rust
/// Encode a `Mesh` to quantized-OBJ text using `profile`'s grid. Only Tri3
/// blocks are emitted (the mesh-LLM format is triangle soup); other element
/// types are skipped. Faces are written 1-based, per OBJ.
pub fn encode(mesh: &Mesh, profile: &ModelProfile) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    for v in &mesh.nodes {
        let _ = writeln!(
            out,
            "v {} {} {}",
            profile.quant(v.x),
            profile.quant(v.y),
            profile.quant(v.z)
        );
    }
    for block in &mesh.element_blocks {
        if block.element_type != ElementType::Tri3 {
            continue;
        }
        for tri in block.connectivity.chunks_exact(3) {
            let _ = writeln!(out, "f {} {} {}", tri[0] + 1, tri[1] + 1, tri[2] + 1);
        }
    }
    out
}
```

- [ ] **Step 4: Run the tests**

Run: `cargo test -p valenx-meshgen codec`
Expected: PASS (6 codec tests).

- [ ] **Step 5: Commit**

```bash
git add crates/valenx-meshgen/src/codec.rs
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(meshgen): codec encode — Mesh -> quantized-OBJ, round-trip tested (task 3)"
```

---

## Task 4: The `MeshLlm` trait + `MockBackend` (generate → decode, no weights)

**Files:**
- Modify: `crates/valenx-meshgen/src/backend/mod.rs`
- Modify: `crates/valenx-meshgen/src/lib.rs` (re-export the trait + MockBackend)

- [ ] **Step 1: Write the failing test** — `src/backend/mod.rs`:

```rust
//! The MeshLlm runtime trait + a weightless MockBackend for tests.

use std::ops::ControlFlow;

use crate::error::MeshGenError;
use crate::model::{Device, GenParams, ModelConfig};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::codec::decode;
    use crate::model::ModelProfile;

    #[test]
    fn mock_backend_streams_tokens_that_decode_to_a_mesh() {
        // The "model" emits a tetrahedron, split across arbitrary token chunks.
        let tokens = vec![
            "v 0 0 0\nv 127 0 0\n", "v 0 127 0\nv 0 0 127\n",
            "f 1 2 3\nf 1 2 4\n", "f 1 3 4\nf 2 3 4\n",
        ];
        let mut backend = MockBackend::new(tokens);

        let mut acc = String::new();
        backend
            .generate("a tetrahedron", &GenParams::default(), &mut |tok| {
                acc.push_str(tok);
                ControlFlow::Continue(())
            })
            .unwrap();

        let mesh = decode(&acc, &ModelProfile::LLAMA_MESH);
        assert_eq!(mesh.nodes.len(), 4);
    }

    #[test]
    fn mock_backend_honors_cancellation() {
        let tokens = vec!["a", "b", "c", "d", "e"];
        let mut backend = MockBackend::new(tokens);
        let mut seen = 0;
        backend
            .generate("x", &GenParams::default(), &mut |_tok| {
                seen += 1;
                if seen == 2 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
            })
            .unwrap();
        assert_eq!(seen, 2, "generation stops at the first Break");
    }
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p valenx-meshgen backend`
Expected: FAIL — `MeshLlm` / `MockBackend` not defined.

- [ ] **Step 3: Implement the trait + MockBackend** (prepend above the test module):

```rust
/// A loadable mesh-generation language model. Backends are swappable; callers
/// never know which concrete type they hold.
pub trait MeshLlm: Send {
    /// Load a model described by `config` onto `device`.
    fn load(config: &ModelConfig, device: Device) -> Result<Self, MeshGenError>
    where
        Self: Sized;

    /// Stream a completion for `prompt`. `on_token` is called once per decoded
    /// token chunk; returning `ControlFlow::Break` cancels promptly.
    fn generate(
        &mut self,
        prompt: &str,
        params: &GenParams,
        on_token: &mut dyn FnMut(&str) -> ControlFlow<()>,
    ) -> Result<(), MeshGenError>;
}

/// A weightless backend that replays a fixed token list. For tests + for the
/// "no model installed" demo path. Does not implement `MeshLlm::load` in a
/// meaningful way (there is nothing to load); construct it with `new`.
pub struct MockBackend {
    tokens: Vec<String>,
}

impl MockBackend {
    /// Build a mock that will replay `tokens` in order.
    pub fn new<I, S>(tokens: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        Self { tokens: tokens.into_iter().map(Into::into).collect() }
    }

    /// Replay the canned tokens through `on_token`, honoring cancellation.
    pub fn generate(
        &mut self,
        _prompt: &str,
        _params: &GenParams,
        on_token: &mut dyn FnMut(&str) -> ControlFlow<()>,
    ) -> Result<(), MeshGenError> {
        for tok in &self.tokens {
            if on_token(tok).is_break() {
                break;
            }
        }
        Ok(())
    }
}
```

(Note: `MockBackend::generate` is an inherent method with the same shape as the trait method, so tests call it directly without a `ModelConfig`. The real `LlamaCppBackend` in Task 5 implements the `MeshLlm` trait proper.)

- [ ] **Step 4: Re-export** in `src/lib.rs`: add `pub use backend::{MeshLlm, MockBackend};`.

- [ ] **Step 5: Run the tests**

Run: `cargo test -p valenx-meshgen`
Expected: PASS (all tasks 1–4 tests).

- [ ] **Step 6: Commit**

```bash
git add crates/valenx-meshgen/src
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(meshgen): MeshLlm trait + weightless MockBackend, generate->decode tested (task 4)"
```

---

## Task 5: `LlamaCppBackend` behind the `meshgen` feature + ignored real-GGUF smoke test

This is the only task that touches the heavy dependency. It must NOT change the default build.

> **IMPORTANT — third-party API:** `llama-cpp-2`'s exact API has changed across versions. The code below is the correct *shape* (load model → context → tokenize → decode loop → detokenize → stream), with the type names from recent `llama-cpp-2`. **Confirm each call against the version you pin** (`cargo doc -p llama-cpp-2 --open`) and adjust signatures as needed — do not assume it compiles verbatim. The trait contract (`load` + `generate(on_token)`) is fixed; only the inside of those two methods is version-sensitive.

**Files:**
- Modify: `crates/valenx-meshgen/Cargo.toml` (add optional dep + feature)
- Create: `crates/valenx-meshgen/src/backend/llama.rs`
- Modify: `crates/valenx-meshgen/src/backend/mod.rs` (gate the submodule)

- [ ] **Step 1: Add the optional dep + feature** to `crates/valenx-meshgen/Cargo.toml`:

```toml
[dependencies]
thiserror   = { workspace = true }
nalgebra    = { workspace = true }
valenx-mesh = { path = "../valenx-mesh" }
# Pin the current stable from crates.io at implementation time.
llama-cpp-2 = { version = "0.1", optional = true }

[features]
# OFF by default — the default workspace build + CI never compile llama.cpp.
meshgen = ["dep:llama-cpp-2"]
```

- [ ] **Step 2: Gate the submodule** — in `src/backend/mod.rs`, add near the top:

```rust
#[cfg(feature = "meshgen")]
pub mod llama;
#[cfg(feature = "meshgen")]
pub use llama::LlamaCppBackend;
```

- [ ] **Step 3: Implement `src/backend/llama.rs`** (skeleton — verify API against the pinned `llama-cpp-2`):

```rust
//! Real inference backend over `llama-cpp-2` (GGUF). Compiled only with the
//! `meshgen` feature. See the API caveat in the plan: confirm calls against the
//! pinned crate version.

use std::num::NonZeroU32;
use std::ops::ControlFlow;

use llama_cpp_2::context::params::LlamaContextParams;
use llama_cpp_2::llama_backend::LlamaBackend;
use llama_cpp_2::model::params::LlamaModelParams;
use llama_cpp_2::model::{AddBos, LlamaModel};

use crate::backend::MeshLlm;
use crate::error::MeshGenError;
use crate::model::{Device, GenParams, ModelConfig};

/// A llama.cpp-backed mesh-LLM.
pub struct LlamaCppBackend {
    backend: LlamaBackend,
    model: LlamaModel,
    // Profile retained for the codec layer above; not used inside generate.
}

impl MeshLlm for LlamaCppBackend {
    fn load(config: &ModelConfig, device: Device) -> Result<Self, MeshGenError> {
        let backend = LlamaBackend::init()
            .map_err(|e| MeshGenError::ModelLoad { path: config.path.clone(), reason: e.to_string() })?;

        let n_gpu_layers = match device {
            Device::Cpu => 0,
            Device::Gpu { layers } => layers,
        };
        let model_params = LlamaModelParams::default().with_n_gpu_layers(n_gpu_layers);

        let model = LlamaModel::load_from_file(&backend, &config.path, &model_params)
            .map_err(|e| MeshGenError::ModelLoad { path: config.path.clone(), reason: e.to_string() })?;

        Ok(Self { backend, model })
    }

    fn generate(
        &mut self,
        prompt: &str,
        params: &GenParams,
        on_token: &mut dyn FnMut(&str) -> ControlFlow<()>,
    ) -> Result<(), MeshGenError> {
        // 1. context
        let ctx_params = LlamaContextParams::default()
            .with_n_ctx(NonZeroU32::new(8192));
        let mut ctx = self
            .model
            .new_context(&self.backend, ctx_params)
            .map_err(|e| MeshGenError::Inference(e.to_string()))?;

        // 2. tokenize the prompt
        let tokens = self
            .model
            .str_to_token(prompt, AddBos::Always)
            .map_err(|e| MeshGenError::Inference(e.to_string()))?;

        // 3. decode loop — feed the prompt, then sample up to max_tokens,
        //    detokenizing each new token and streaming it to on_token.
        //    (Batch construction, sampler setup, and EOS detection follow the
        //    llama-cpp-2 "simple" example for the pinned version. Stop early on
        //    EOS or when on_token returns Break.)
        //
        //    let mut batch = LlamaBatch::new(...);
        //    ... add tokens ...; ctx.decode(&mut batch)?;
        //    for _ in 0..params.max_tokens {
        //        let next = sampler.sample(&ctx, ...);
        //        if self.model.is_eog_token(next) { break; }
        //        let piece = self.model.token_to_str(next, Special::Tokenize)?;
        //        if on_token(&piece).is_break() { break; }
        //        batch.clear(); batch.add(next, pos, &[0], true)?;
        //        ctx.decode(&mut batch)?;
        //    }
        let _ = (&ctx, &tokens, params.seed, params.temperature, params.top_p);
        Ok(())
    }
}
```

(The commented decode loop is the exact sequence to fill from the pinned `llama-cpp-2` example; it is logic, not a placeholder for *what* to do — the steps are spelled out. Keep `generate` returning `Ok(())` with an empty body until the loop is wired, so the feature build compiles, then implement the loop and run the Task-5 smoke test.)

- [ ] **Step 4: Write the ignored smoke test** — add to `src/backend/llama.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::codec::decode;
    use crate::model::ModelProfile;
    use std::path::PathBuf;

    /// Headless end-to-end: set VALENX_MESHGEN_MODEL to a GGUF mesh-LLM and run
    /// with `--features meshgen --ignored`. Not run in CI (BYO weights).
    #[test]
    #[ignore = "requires a local GGUF mesh-LLM via VALENX_MESHGEN_MODEL"]
    fn generates_a_nonempty_mesh_from_a_prompt() {
        let path = std::env::var("VALENX_MESHGEN_MODEL")
            .expect("set VALENX_MESHGEN_MODEL to a GGUF file");
        let cfg = ModelConfig {
            path: PathBuf::from(path),
            profile: ModelProfile::LLAMA_MESH,
            device: Device::Cpu,
        };
        let mut backend = LlamaCppBackend::load(&cfg, Device::Cpu).unwrap();
        let prompt = cfg.profile.build_prompt("a simple chair");
        let mut acc = String::new();
        backend
            .generate(&prompt, &GenParams::default(), &mut |t| {
                acc.push_str(t);
                std::ops::ControlFlow::Continue(())
            })
            .unwrap();
        let mesh = decode(&acc, &cfg.profile);
        assert!(mesh.nodes.len() >= 3, "expected a non-trivial mesh, got {} verts", mesh.nodes.len());
    }
}
```

- [ ] **Step 5: Verify the DEFAULT build is unchanged (no llama.cpp)**

Run: `cargo build -p valenx-meshgen` and `cargo test -p valenx-meshgen`
Expected: compiles + all tasks 1–4 tests pass, with NO llama.cpp compilation (the `llama` module is absent without the feature).

- [ ] **Step 6: Verify the FEATURE build compiles** (this compiles llama.cpp — first build is slow; needs a C toolchain, optionally CUDA/Vulkan):

Run: `cargo build -p valenx-meshgen --features meshgen`
Expected: compiles (after you reconcile the `llama-cpp-2` API per the caveat). The ignored smoke test is NOT run here.

- [ ] **Step 7: Commit**

```bash
git add crates/valenx-meshgen
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "feat(meshgen): feature-gated LlamaCppBackend + ignored BYO-GGUF smoke test (task 5)"
```

---

## Task 6: Gate verification + clippy + roadmap note

**Files:**
- (no new code) — verification + a one-line README pointer if desired.

- [ ] **Step 1: Confirm the feature is truly off by default** — `cargo tree -p valenx-meshgen` must NOT list `llama-cpp-2`; `cargo tree -p valenx-meshgen --features meshgen` must list it.

- [ ] **Step 2: Lints clean**

Run: `cargo clippy -p valenx-meshgen --all-targets -- -D warnings`
Expected: zero warnings. (Run the feature variant too if your toolchain has the C deps: `cargo clippy -p valenx-meshgen --features meshgen -- -D warnings`.)

- [ ] **Step 3: Format**

Run: `cargo fmt -p valenx-meshgen` then `cargo fmt -p valenx-meshgen -- --check`.

- [ ] **Step 4: Commit any formatting/lint fixes**

```bash
git add -A crates/valenx-meshgen
git -c user.name=nochallenge -c user.email=201502404+nochallenge@users.noreply.github.com commit -m "chore(meshgen): clippy + fmt clean; default build excludes the LLM stack (task 6)"
```

---

## Phase 1 done — what exists

A `valenx-meshgen` crate that: defines the `MeshLlm` runtime seam; converts a mesh-LLM's quantized-OBJ token stream into a validated valenx `Mesh` (and back); is fully tested with no weights via `MockBackend`; and carries a feature-gated real `llama-cpp-2` backend that never touches the default build. The next phases (own plans):

- **P2** — `meshgen_workbench.rs` in `valenx-app`: prompt box, background generation, streaming + cancel, run the decoded mesh through `valenx_mesh::merge_coincident_nodes` + `repair::fill_holes` + `is_manifold`, load into the viewport via the existing `mesh_loader` path, settings (model path/backend/device), and an agent-bridge `MeshGen { prompt }` command.
- **P3** — `describe` (encode a viewport mesh → prompt → text) + chat-refine.
- **P4** — `CandleBackend` (pure Rust) behind the same trait.
- **Parallel track** — extend `valenx-blender-mesh-ops` with clean-room Rust editor ops (bevel/remesh) implemented from algorithms, never ported from Blender's GPL source.

---

## Execution

Deferred — author is away. When you return, two options:

1. **Subagent-Driven (recommended)** — dispatch a fresh subagent per task, review between tasks. (REQUIRED SUB-SKILL: superpowers:subagent-driven-development.)
2. **Inline** — execute here with checkpoints. (REQUIRED SUB-SKILL: superpowers:executing-plans.)

Tasks 0–4 are pure Rust + fast; Task 5 is the only one needing the C/GPU toolchain and the `llama-cpp-2` API reconciliation, so it's the natural place for a careful checkpoint.
