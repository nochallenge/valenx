//! Frame metadata + aggregate material-removal report types.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Per-frame snapshot.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FrameMetadata {
    /// Frame index (0-based).
    pub index: usize,
    /// Symbolic time `t ∈ [0, 1]`.
    pub t: f64,
    /// Tool position at this frame (mm, stock-local).
    pub position: Vector3<f64>,
    /// Total cut volume since `t = 0` (mm³).
    pub cut_volume_mm3: f64,
    /// Swarf volume — equal to `cut_volume_mm3` for solid stock (no
    /// pre-existing holes).
    pub swarf_volume_mm3: f64,
    /// Material removal rate over the previous frame interval
    /// (mm³ / unit-t). v1 has no notion of real machining time so the
    /// unit is `Δt = 1/(n_frames - 1)`.
    pub mrr_mm3_per_unit_t: f64,
}

/// Aggregate report across every frame.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MaterialRemovalReport {
    /// Frame-by-frame metadata.
    pub frames: Vec<FrameMetadata>,
    /// Total cut volume = `frames.last().cut_volume_mm3`.
    pub total_cut_volume_mm3: f64,
    /// Tool deflection (mm). v1 always 0.0 — placeholder until the
    /// stiffness model lands.
    pub tool_deflection_mm: f64,
}

impl MaterialRemovalReport {
    /// Convenience accessor — peak per-frame MRR.
    pub fn peak_mrr(&self) -> f64 {
        self.frames
            .iter()
            .map(|f| f.mrr_mm3_per_unit_t)
            .fold(0.0_f64, f64::max)
    }
}
