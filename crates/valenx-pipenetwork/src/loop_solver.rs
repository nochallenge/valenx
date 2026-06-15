//! Loops and the per-loop Hardy-Cross correction.
//!
//! A [`Loop`] is an ordered, oriented set of pipes forming a closed circuit
//! in the network graph. Each member carries an [`Orientation`] recording
//! whether the pipe's own reference direction agrees (`Forward`) or
//! disagrees (`Reverse`) with the direction in which the loop is traversed.
//!
//! For a loop `L`, define for each member the loop-aligned flow
//! `s_i * q_i`, where `s_i = +1` for `Forward` and `s_i = -1` for
//! `Reverse`. The Hardy-Cross loop correction is
//!
//! ```text
//! dQ = - sum_i ( k_i * (s_i q_i) * |s_i q_i| )  /  sum_i ( 2 k_i * |s_i q_i| )
//!    = - sum_i ( s_i * k_i * q_i * |q_i| )       /  sum_i ( 2 k_i * |q_i| )
//! ```
//!
//! (using `|s_i q_i| = |q_i|` since `|s_i| = 1`). The correction `dQ` is
//! then applied to every member as `q_i <- q_i + s_i * dQ`, which keeps the
//! flow continuous around the loop. At a balanced loop the numerator — the
//! signed sum of head losses around the loop — is zero, so `dQ` is zero and
//! the flows stop changing.

use crate::error::NetworkError;
use crate::pipe::Pipe;
use serde::{Deserialize, Serialize};

/// Whether a pipe's reference direction agrees with the loop traversal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Orientation {
    /// The pipe's positive-`q` direction matches the loop's traversal
    /// direction; its loop sign `s` is `+1`.
    Forward,
    /// The pipe's positive-`q` direction is opposite to the loop's
    /// traversal direction; its loop sign `s` is `-1`.
    Reverse,
}

impl Orientation {
    /// The loop sign `s` for this orientation: `+1.0` for [`Forward`],
    /// `-1.0` for [`Reverse`].
    ///
    /// [`Forward`]: Orientation::Forward
    /// [`Reverse`]: Orientation::Reverse
    pub fn sign(self) -> f64 {
        match self {
            Orientation::Forward => 1.0,
            Orientation::Reverse => -1.0,
        }
    }
}

/// One membership of a pipe in a loop: which pipe, and how it is oriented
/// relative to the loop's traversal direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct LoopMember {
    /// Index of the pipe within the owning network's pipe list.
    pub pipe: usize,
    /// Orientation of the pipe relative to loop traversal.
    pub orientation: Orientation,
}

impl LoopMember {
    /// A forward-oriented membership of the pipe at `index`.
    pub fn forward(index: usize) -> Self {
        Self {
            pipe: index,
            orientation: Orientation::Forward,
        }
    }

    /// A reverse-oriented membership of the pipe at `index`.
    pub fn reverse(index: usize) -> Self {
        Self {
            pipe: index,
            orientation: Orientation::Reverse,
        }
    }
}

/// A closed loop in the network: a named, ordered list of oriented pipe
/// memberships.
///
/// The order of members is irrelevant to the Hardy-Cross correction (the
/// formula is a sum), but is preserved so callers can reason about the
/// physical circuit. Loops are validated against the network at solve time.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Loop {
    /// A human-readable label for diagnostics (e.g. `"loop-1"`).
    pub name: String,
    /// The oriented pipe memberships that make up the loop.
    pub members: Vec<LoopMember>,
}

impl Loop {
    /// Build a loop from a name and its members.
    ///
    /// # Errors
    ///
    /// Returns [`NetworkError::EmptyLoop`] if `members` is empty.
    pub fn new(name: impl Into<String>, members: Vec<LoopMember>) -> Result<Self, NetworkError> {
        let name = name.into();
        if members.is_empty() {
            return Err(NetworkError::EmptyLoop(name));
        }
        Ok(Self { name, members })
    }

    /// Validate every member index against a network of `pipe_count` pipes.
    ///
    /// # Errors
    ///
    /// Returns [`NetworkError::UnknownPipe`] for the first member whose
    /// pipe index is `>= pipe_count`.
    pub fn validate_against(&self, pipe_count: usize) -> Result<(), NetworkError> {
        for member in &self.members {
            if member.pipe >= pipe_count {
                return Err(NetworkError::UnknownPipe {
                    index: member.pipe,
                    count: pipe_count,
                });
            }
        }
        Ok(())
    }

    /// Signed sum of head losses around this loop, `sum_i s_i * k_i q_i|q_i|`.
    ///
    /// This is the numerator of the Hardy-Cross correction. It is the
    /// quantity driven to zero at convergence: a balanced loop has no net
    /// head loss around its circuit.
    ///
    /// `pipes` must be the owning network's pipe slice; member indices are
    /// assumed valid (call [`Loop::validate_against`] first).
    pub fn head_loss_sum(&self, pipes: &[Pipe]) -> f64 {
        self.members
            .iter()
            .map(|m| m.orientation.sign() * pipes[m.pipe].head_loss())
            .sum()
    }

    /// Sum of head-loss slopes around this loop, `sum_i 2 k_i |q_i|`.
    ///
    /// This is the denominator of the Hardy-Cross correction. The loop sign
    /// drops out because the slope depends only on `|q_i|`.
    ///
    /// `pipes` must be the owning network's pipe slice; member indices are
    /// assumed valid.
    pub fn slope_sum(&self, pipes: &[Pipe]) -> f64 {
        self.members
            .iter()
            .map(|m| pipes[m.pipe].head_loss_slope())
            .sum()
    }

    /// Compute the Hardy-Cross loop correction `dQ` for the current flows.
    ///
    /// Returns
    ///
    /// ```text
    /// dQ = - head_loss_sum / slope_sum
    /// ```
    ///
    /// or `None` when `slope_sum` is zero — which happens only when every
    /// pipe in the loop is simultaneously loss-free (`k = 0`) or stationary
    /// (`q = 0`), leaving the correction undefined. Callers treat that as a
    /// degenerate loop.
    ///
    /// `pipes` must be the owning network's pipe slice; member indices are
    /// assumed valid.
    pub fn correction(&self, pipes: &[Pipe]) -> Option<f64> {
        let denominator = self.slope_sum(pipes);
        if denominator == 0.0 {
            return None;
        }
        Some(-self.head_loss_sum(pipes) / denominator)
    }

    /// Apply a loop correction `dq` to every member's flow in `pipes`,
    /// honouring each member's orientation: `q_i <- q_i + s_i * dq`.
    ///
    /// `pipes` must be the owning network's pipe slice; member indices are
    /// assumed valid.
    pub fn apply_correction(&self, pipes: &mut [Pipe], dq: f64) {
        for member in &self.members {
            pipes[member.pipe].q += member.orientation.sign() * dq;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-12;

    /// A parallel two-pipe loop: pipe 0 traversed forward, pipe 1 reverse.
    /// k0=1, k1=4; initial flows q0=q1=1.5.
    fn parallel_loop() -> (Vec<Pipe>, Loop) {
        let pipes = vec![Pipe::new(1.0, 1.5).unwrap(), Pipe::new(4.0, 1.5).unwrap()];
        let lp = Loop::new(
            "loop-1",
            vec![LoopMember::forward(0), LoopMember::reverse(1)],
        )
        .unwrap();
        (pipes, lp)
    }

    #[test]
    fn orientation_signs() {
        assert!((Orientation::Forward.sign() - 1.0).abs() < EPS);
        assert!((Orientation::Reverse.sign() + 1.0).abs() < EPS);
    }

    #[test]
    fn head_loss_sum_uses_orientation() {
        let (pipes, lp) = parallel_loop();
        // s0 k0 q0|q0| + s1 k1 q1|q1|
        //   = (+1)(1)(1.5)(1.5) + (-1)(4)(1.5)(1.5) = 2.25 - 9 = -6.75
        assert!((lp.head_loss_sum(&pipes) - (-6.75)).abs() < EPS);
    }

    #[test]
    fn slope_sum_ignores_orientation() {
        let (pipes, lp) = parallel_loop();
        // 2 k0 |q0| + 2 k1 |q1| = 2*1*1.5 + 2*4*1.5 = 3 + 12 = 15
        assert!((lp.slope_sum(&pipes) - 15.0).abs() < EPS);
    }

    #[test]
    fn correction_matches_hand_computation() {
        let (pipes, lp) = parallel_loop();
        // dQ = -(-6.75)/15 = 0.45
        let dq = lp.correction(&pipes).unwrap();
        assert!((dq - 0.45).abs() < EPS);
    }

    #[test]
    fn correction_equals_negative_loss_over_slope() {
        // Independent of the specific numbers: dQ = -head_loss_sum/slope_sum.
        let pipes = vec![Pipe::new(2.0, 0.8).unwrap(), Pipe::new(0.5, 1.3).unwrap()];
        let lp = Loop::new("l", vec![LoopMember::forward(0), LoopMember::reverse(1)]).unwrap();
        let dq = lp.correction(&pipes).unwrap();
        let expected = -lp.head_loss_sum(&pipes) / lp.slope_sum(&pipes);
        assert!((dq - expected).abs() < EPS);
    }

    #[test]
    fn apply_correction_respects_orientation() {
        let (mut pipes, lp) = parallel_loop();
        lp.apply_correction(&mut pipes, 0.45);
        // forward member gains +dQ, reverse member gains -dQ.
        assert!((pipes[0].q - (1.5 + 0.45)).abs() < EPS);
        assert!((pipes[1].q - (1.5 - 0.45)).abs() < EPS);
    }

    #[test]
    fn apply_correction_preserves_flow_split_sum() {
        // q0 + q1 is invariant under the loop correction (continuity).
        let (mut pipes, lp) = parallel_loop();
        let before = pipes[0].q + pipes[1].q;
        lp.apply_correction(&mut pipes, 0.37);
        let after = pipes[0].q + pipes[1].q;
        assert!((before - after).abs() < EPS);
    }

    #[test]
    fn degenerate_loop_correction_is_none() {
        // All-zero flow => slope sum 0 => correction undefined.
        let pipes = vec![Pipe::new(1.0, 0.0).unwrap(), Pipe::new(4.0, 0.0).unwrap()];
        let lp = Loop::new("l", vec![LoopMember::forward(0), LoopMember::reverse(1)]).unwrap();
        assert!(lp.correction(&pipes).is_none());
    }

    #[test]
    fn empty_loop_is_rejected() {
        let err = Loop::new("oops", Vec::new()).unwrap_err();
        assert!(matches!(err, NetworkError::EmptyLoop(name) if name == "oops"));
    }

    #[test]
    fn validate_against_flags_out_of_range_member() {
        let lp = Loop::new("l", vec![LoopMember::forward(0), LoopMember::forward(5)]).unwrap();
        let err = lp.validate_against(2).unwrap_err();
        assert!(matches!(
            err,
            NetworkError::UnknownPipe { index: 5, count: 2 }
        ));
        // In-range is fine.
        assert!(lp.validate_against(6).is_ok());
    }
}
