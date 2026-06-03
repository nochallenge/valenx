//! UI panel state envelope for the CAMotics animation.

use crate::animation::Animation;
use crate::error::CamoticsError;
use crate::report::MaterialRemovalReport;

/// Workbench-panel state.
pub struct CamoticsPanelState {
    /// Optional animation — `None` until the user loads stock + toolpath.
    pub animation: Option<Animation>,
    /// Current playback `t ∈ [0, 1]`.
    pub current_t: f64,
    /// Cached report after a Run.
    pub last_report: Option<MaterialRemovalReport>,
    /// Last status message.
    pub last_status: Option<String>,
    /// Last error message.
    pub last_error: Option<String>,
}

impl Default for CamoticsPanelState {
    fn default() -> Self {
        Self {
            animation: None,
            current_t: 0.0,
            last_report: None,
            last_status: None,
            last_error: None,
        }
    }
}

impl CamoticsPanelState {
    /// Empty state.
    pub fn new() -> Self {
        Self::default()
    }

    /// Install a new animation; resets `current_t` and clears cache.
    pub fn set_animation(&mut self, anim: Animation) {
        self.current_t = 0.0;
        self.last_report = None;
        self.animation = Some(anim);
        self.last_status = Some("animation loaded".into());
        self.last_error = None;
    }

    /// Move playback. Clamps to `[0, 1]`.
    pub fn seek(&mut self, t: f64) {
        self.current_t = t.clamp(0.0, 1.0);
    }

    /// Run the full report — replays every frame and stashes it.
    pub fn compute_report(&mut self) -> Result<(), CamoticsError> {
        let anim = self
            .animation
            .as_ref()
            .ok_or_else(|| CamoticsError::BadParameter {
                name: "animation",
                reason: "no animation loaded".into(),
            })?;
        let r = anim.material_removal_report()?;
        self.last_status = Some(format!(
            "{} frames, total cut {:.2} mm³ (peak MRR {:.2})",
            r.frames.len(),
            r.total_cut_volume_mm3,
            r.peak_mrr()
        ));
        self.last_error = None;
        self.last_report = Some(r);
        Ok(())
    }

    /// Record a status message.
    pub fn set_status(&mut self, msg: impl Into<String>) {
        self.last_status = Some(msg.into());
        self.last_error = None;
    }

    /// Record an error message.
    pub fn set_error(&mut self, msg: impl Into<String>) {
        self.last_error = Some(msg.into());
        self.last_status = None;
    }
}
