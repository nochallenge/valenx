//! Model configuration + the quantization/prompt profile.

use std::path::PathBuf;

/// Where inference runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Device {
    /// CPU only.
    Cpu,
    /// Offload `layers` transformer layers to the GPU (llama.cpp `n_gpu_layers`).
    Gpu {
        /// Number of transformer layers to offload to the GPU.
        layers: u32,
    },
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
        Self {
            max_tokens: 8192,
            temperature: 0.7,
            top_p: 0.9,
            seed: 0,
        }
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
