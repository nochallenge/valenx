//! Per-frame atomic coordinates from MD output.
//!
//! A `Trajectory` is the canonical Valenx representation of an MD
//! simulation snapshot stream: an outer per-frame `Vec` of inner
//! per-atom `Vec<Vector3<f64>>`. The atom count is constant across
//! frames; the [`Trajectory::validate`] helper enforces that contract
//! after a multi-frame load (e.g. from the [`crate::format::dcd`]
//! reader) so downstream code can rely on it.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Per-frame atomic coordinate stack from an MD trajectory.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct Trajectory {
    /// Source identifier (e.g. DCD file basename, simulation name).
    pub id: String,
    /// Per-frame `Vec<Vector3<f64>>` — outer length is frame count,
    /// inner length is atom count (constant across frames).
    pub frames: Vec<Vec<Vector3<f64>>>,
}

/// Errors surfaced by [`Trajectory::validate`].
#[derive(Debug, Error, PartialEq, Eq)]
pub enum TrajectoryError {
    /// At least one frame's atom count differs from the first frame's.
    /// `expected` is the atom count anchored on frame 0; `frame` is
    /// the 0-indexed offending frame; `found` is its atom count.
    #[error(
        "trajectory frames have inconsistent atom counts: \
         frame 0 has {expected} atoms, frame {frame} has {found}"
    )]
    Inconsistent {
        /// Atom count seen in frame 0.
        expected: usize,
        /// 0-indexed offending frame.
        frame: usize,
        /// Atom count seen in the offending frame.
        found: usize,
    },
}

impl Trajectory {
    /// Build a trajectory from an id and an owned frame stack. No
    /// validation runs here — call [`Trajectory::validate`] after
    /// construction if you need the consistent-atom-count contract.
    pub fn new(id: impl Into<String>, frames: Vec<Vec<Vector3<f64>>>) -> Self {
        Self {
            id: id.into(),
            frames,
        }
    }

    /// Number of frames in the trajectory.
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Atom count of frame 0, or `None` if the trajectory has no frames.
    pub fn atom_count(&self) -> Option<usize> {
        self.frames.first().map(|f| f.len())
    }

    /// Bounds-checked frame access. Returns `None` if `i` is out of
    /// range; otherwise a borrowed slice of the frame's atom
    /// coordinates.
    pub fn frame(&self, i: usize) -> Option<&[Vector3<f64>]> {
        self.frames.get(i).map(|v| v.as_slice())
    }

    /// Validate that every frame has the same atom count as frame 0.
    /// Empty trajectories validate trivially. Returns
    /// [`TrajectoryError::Inconsistent`] on the first mismatch.
    pub fn validate(&self) -> Result<(), TrajectoryError> {
        let mut frames = self.frames.iter().enumerate();
        let Some((_, first)) = frames.next() else {
            return Ok(());
        };
        let expected = first.len();
        for (idx, frame) in frames {
            if frame.len() != expected {
                return Err(TrajectoryError::Inconsistent {
                    expected,
                    frame: idx,
                    found: frame.len(),
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_trajectory_has_zero_frames() {
        let t = Trajectory::default();
        assert_eq!(t.frame_count(), 0);
        assert_eq!(t.atom_count(), None);
    }

    #[test]
    fn frame_returns_some_for_valid_index() {
        let frames = vec![
            vec![Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 5.0, 6.0)],
            vec![Vector3::new(7.0, 8.0, 9.0), Vector3::new(0.0, 0.0, 0.0)],
        ];
        let t = Trajectory::new("traj", frames);
        let frame0 = t.frame(0).expect("frame 0 in range");
        assert_eq!(frame0.len(), 2);
        assert_eq!(frame0[0], Vector3::new(1.0, 2.0, 3.0));
        // Out-of-range returns None rather than panicking.
        assert!(t.frame(5).is_none());
    }

    #[test]
    fn validate_passes_consistent_frames() {
        // Three frames, each carrying the same 3-atom coordinate
        // count. Validation must succeed.
        let frames = vec![
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
            ],
            vec![
                Vector3::new(0.1, 0.0, 0.0),
                Vector3::new(1.1, 0.0, 0.0),
                Vector3::new(0.1, 1.0, 0.0),
            ],
            vec![
                Vector3::new(0.2, 0.0, 0.0),
                Vector3::new(1.2, 0.0, 0.0),
                Vector3::new(0.2, 1.0, 0.0),
            ],
        ];
        let t = Trajectory::new("ok", frames);
        assert!(t.validate().is_ok());
    }

    #[test]
    fn validate_rejects_mismatched_frame_count() {
        // Frame 1 carries 2 atoms instead of 3 — the validator must
        // surface the offending index + counts.
        let frames = vec![
            vec![
                Vector3::new(0.0, 0.0, 0.0),
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.0, 1.0, 0.0),
            ],
            vec![Vector3::new(0.1, 0.0, 0.0), Vector3::new(1.1, 0.0, 0.0)],
        ];
        let t = Trajectory::new("bad", frames);
        let err = t.validate().expect_err("inconsistent frame count");
        assert_eq!(
            err,
            TrajectoryError::Inconsistent {
                expected: 3,
                frame: 1,
                found: 2,
            }
        );
    }
}
