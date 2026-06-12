//! Keyframe animation data model + sampler.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::error::AnimateError;
use crate::tween::TweenMode;

/// Joint id type — matches `valenx_assembly::Joint::id` (`usize`).
pub type JointId = usize;

/// One keyframe: a snapshot of joint parameters at a given time, plus
/// the tween mode used to ease *into* this keyframe from its
/// predecessor.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Keyframe {
    /// Seconds from the start of the animation.
    pub time: f64,
    /// Per-joint parameter target. Joints missing from this map keep
    /// whatever value they had at the previous keyframe (i.e. sparse
    /// keyframes are allowed).
    pub joint_parameters: HashMap<JointId, f64>,
    /// Tween mode applied when interpolating from the previous
    /// keyframe to this one. The first keyframe's `tween` is
    /// ignored.
    pub tween: TweenMode,
}

impl Keyframe {
    /// Build a keyframe at `time` with linear tween and no joint
    /// parameters.
    pub fn at(time: f64) -> Self {
        Self {
            time,
            joint_parameters: HashMap::new(),
            tween: TweenMode::Linear,
        }
    }

    /// Convenience: set a joint parameter (chainable).
    pub fn with_joint(mut self, id: JointId, value: f64) -> Self {
        self.joint_parameters.insert(id, value);
        self
    }

    /// Convenience: set tween mode (chainable).
    pub fn with_tween(mut self, tween: TweenMode) -> Self {
        self.tween = tween;
        self
    }
}

/// One animation = an ordered list of keyframes.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Animation {
    /// Optional user-facing name.
    pub name: String,
    /// Keyframes in monotonically-increasing time order.
    pub keyframes: Vec<Keyframe>,
}

impl Animation {
    /// Empty animation.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a keyframe; returns `Err(NotMonotonic)` if it would
    /// break the monotonic-time invariant.
    pub fn push(&mut self, kf: Keyframe) -> Result<(), AnimateError> {
        if let Some(last) = self.keyframes.last() {
            if kf.time < last.time {
                return Err(AnimateError::NotMonotonic {
                    reason: format!("incoming t={} < last t={}", kf.time, last.time),
                });
            }
        }
        self.keyframes.push(kf);
        Ok(())
    }

    /// Total duration = last keyframe time. 0.0 when empty.
    pub fn duration(&self) -> f64 {
        self.keyframes.last().map(|k| k.time).unwrap_or(0.0)
    }

    /// Sample the joint parameter map at time `t`.
    ///
    /// Each joint id ever mentioned in the keyframes gets a value;
    /// the value is the linear-or-eased interpolation between the two
    /// bracketing keyframes that mention that joint. Joints whose
    /// last-mention is before `t` carry their final value forward.
    /// Joints whose first-mention is after `t` are absent from the
    /// output.
    pub fn sample(&self, t: f64) -> Result<HashMap<JointId, f64>, AnimateError> {
        if self.keyframes.is_empty() {
            return Err(AnimateError::Empty);
        }
        // Collect every joint id ever mentioned.
        let mut all_ids: std::collections::BTreeSet<JointId> = std::collections::BTreeSet::new();
        for kf in &self.keyframes {
            for k in kf.joint_parameters.keys() {
                all_ids.insert(*k);
            }
        }
        let mut out = HashMap::new();
        for id in all_ids {
            // Find the (k_prev, k_next) bracketing this joint at time t.
            let mut prev: Option<&Keyframe> = None;
            let mut next: Option<&Keyframe> = None;
            for kf in &self.keyframes {
                if !kf.joint_parameters.contains_key(&id) {
                    continue;
                }
                if kf.time <= t {
                    prev = Some(kf);
                } else if next.is_none() {
                    next = Some(kf);
                    break;
                }
            }
            let v = match (prev, next) {
                (Some(p), Some(n)) => {
                    let span = n.time - p.time;
                    let local = if span.abs() < f64::EPSILON {
                        0.0
                    } else {
                        ((t - p.time) / span).clamp(0.0, 1.0)
                    };
                    let eased = n.tween.apply(local);
                    let pv = *p.joint_parameters.get(&id).unwrap();
                    let nv = *n.joint_parameters.get(&id).unwrap();
                    pv + (nv - pv) * eased
                }
                (Some(p), None) => *p.joint_parameters.get(&id).unwrap(),
                (None, _) => continue, // joint hasn't appeared yet
            };
            out.insert(id, v);
        }
        Ok(out)
    }
}

/// One frame in a playback: the time the frame represents and the
/// resolved joint-parameter map.
#[derive(Clone, Debug, PartialEq)]
pub struct FrameSnapshot {
    /// Seconds from animation start.
    pub time: f64,
    /// Resolved joint id → parameter.
    pub joint_parameters: HashMap<JointId, f64>,
}

impl Animation {
    /// Drive the animation over `[start_t, end_t]` at `fps` and
    /// return the per-frame joint-parameter snapshots.
    ///
    /// Does **not** actually mutate or render the assembly — that's
    /// left to the caller (the desktop shell pipes each snapshot
    /// into `valenx_assembly::kinematics::apply_all_joints` against
    /// a clone of the assembly).
    pub fn frames(
        &self,
        start_t: f64,
        end_t: f64,
        fps: u32,
    ) -> Result<Vec<FrameSnapshot>, AnimateError> {
        if fps == 0 {
            return Err(AnimateError::BadParameter {
                name: "fps",
                reason: "must be > 0".into(),
            });
        }
        if end_t <= start_t {
            return Err(AnimateError::BadParameter {
                name: "end_t",
                reason: "must be > start_t".into(),
            });
        }
        let span_frames = (end_t - start_t) * fps as f64;
        // Bound the frame count so a huge span/fps can't allocate gigabytes or
        // hang the sample loop (~4.6 h at 60 fps; mirrors valenx-optimize's
        // size caps). A non-finite product (e.g. end_t = f64::MAX) is rejected too.
        const MAX_FRAMES: f64 = 1_000_000.0;
        if !span_frames.is_finite() || span_frames > MAX_FRAMES {
            return Err(AnimateError::BadParameter {
                name: "end_t",
                reason: format!(
                    "(end_t - start_t) * fps = {span_frames} exceeds the {MAX_FRAMES} frame cap"
                ),
            });
        }
        let count = span_frames.ceil() as usize;
        let mut out = Vec::with_capacity(count + 1);
        for i in 0..=count {
            let t = start_t + (i as f64) / (fps as f64);
            let t = t.min(end_t);
            let map = self.sample(t)?;
            out.push(FrameSnapshot {
                time: t,
                joint_parameters: map,
            });
        }
        Ok(out)
    }

    /// Apply one frame's joint parameters to an assembly and run
    /// `apply_all_joints` to update part transforms. Returns the
    /// mutated assembly clone (the input `base` is left unmodified).
    pub fn apply_frame(
        base: &valenx_assembly::Assembly,
        frame: &FrameSnapshot,
    ) -> Result<valenx_assembly::Assembly, AnimateError> {
        let mut a = base.clone();
        // Overwrite each named joint's parameter.
        for j in a.joints.iter_mut() {
            if let Some(v) = frame.joint_parameters.get(&j.id) {
                j.parameter = *v;
            }
        }
        valenx_assembly::kinematics::apply_all_joints(&mut a)?;
        Ok(a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_sample_errors() {
        let a = Animation::new();
        assert!(a.sample(0.0).is_err());
    }

    #[test]
    fn linear_interpolation_midpoint() {
        let mut a = Animation::new();
        a.push(Keyframe::at(0.0).with_joint(0, 0.0)).unwrap();
        a.push(Keyframe::at(2.0).with_joint(0, 10.0)).unwrap();
        let s = a.sample(1.0).unwrap();
        assert!((s.get(&0).copied().unwrap() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn nonmonotonic_errors() {
        let mut a = Animation::new();
        a.push(Keyframe::at(1.0)).unwrap();
        assert!(a.push(Keyframe::at(0.5)).is_err());
    }

    #[test]
    fn frames_count_matches_fps() {
        let mut a = Animation::new();
        a.push(Keyframe::at(0.0).with_joint(0, 0.0)).unwrap();
        a.push(Keyframe::at(1.0).with_joint(0, 1.0)).unwrap();
        let frames = a.frames(0.0, 1.0, 10).unwrap();
        // 10 fps over 1 s with inclusive endpoints → 11 frames
        assert_eq!(frames.len(), 11);
        assert!((frames.last().unwrap().joint_parameters[&0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn fps_zero_errors() {
        let mut a = Animation::new();
        a.push(Keyframe::at(0.0).with_joint(0, 0.0)).unwrap();
        a.push(Keyframe::at(1.0).with_joint(0, 1.0)).unwrap();
        assert!(a.frames(0.0, 1.0, 0).is_err());
    }

    #[test]
    fn frames_rejects_huge_span_not_oom() {
        // span × fps = 10^10 would otherwise allocate ~240 GB and hang the
        // sample loop; it must be rejected.
        let mut a = Animation::new();
        a.push(Keyframe::at(0.0).with_joint(0, 0.0)).unwrap();
        a.push(Keyframe::at(1.0).with_joint(0, 1.0)).unwrap();
        assert!(
            a.frames(0.0, 1e7, 1000).is_err(),
            "huge span must be rejected"
        );
        // A non-finite product (end_t = f64::MAX) is rejected too.
        assert!(a.frames(0.0, f64::MAX, 1000).is_err());
    }

    #[test]
    fn carry_forward_after_last_keyframe() {
        let mut a = Animation::new();
        a.push(Keyframe::at(0.0).with_joint(0, 1.0)).unwrap();
        a.push(Keyframe::at(1.0).with_joint(0, 5.0)).unwrap();
        let s = a.sample(2.0).unwrap();
        assert!((s.get(&0).copied().unwrap() - 5.0).abs() < 1e-9);
    }

    #[test]
    fn easing_in_out_curves() {
        let mut a = Animation::new();
        a.push(Keyframe::at(0.0).with_joint(0, 0.0)).unwrap();
        a.push(
            Keyframe::at(1.0)
                .with_joint(0, 10.0)
                .with_tween(TweenMode::EaseInOut),
        )
        .unwrap();
        let s = a.sample(0.5).unwrap();
        // Smooth-step at 0.5 = 0.5 → midpoint same as linear
        assert!((s.get(&0).copied().unwrap() - 5.0).abs() < 1e-6);
    }
}
