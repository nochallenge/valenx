//! Solved quantities: per-member axial forces and per-support reactions.
//!
//! [`crate::solver::solve`] returns a [`TrussSolution`] bundling the
//! solved axial force in every member and the reaction at every support.
//! The convenience methods classify each member as tension / compression
//! and let a caller (or the test-suite) re-check global equilibrium
//! against the applied loads.

use serde::{Deserialize, Serialize};

use crate::model::Truss;

/// Classification of a member's axial state from the sign of its force.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AxialState {
    /// Positive axial force: the member is stretched (pulls its joints
    /// together).
    Tension,
    /// Negative axial force: the member is squashed (pushes its joints
    /// apart).
    Compression,
    /// Force is (within tolerance) zero — a *zero-force member*.
    Zero,
}

/// The solved axial force in one member.
///
/// `force` is **tension positive**: `> 0` tension, `< 0` compression.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MemberForce {
    /// Index of this member in [`Truss::members`].
    pub member: usize,
    /// First end-node index.
    pub a: usize,
    /// Second end-node index.
    pub b: usize,
    /// Axial force, tension positive (same unit as the applied loads).
    pub force: f64,
    /// Member length (same unit as the node coordinates).
    pub length: f64,
}

impl MemberForce {
    /// Classify the member as [`AxialState::Tension`],
    /// [`AxialState::Compression`], or [`AxialState::Zero`], treating any
    /// `|force| ≤ tol` as zero.
    #[must_use]
    pub fn state_with_tol(&self, tol: f64) -> AxialState {
        if self.force.abs() <= tol {
            AxialState::Zero
        } else if self.force > 0.0 {
            AxialState::Tension
        } else {
            AxialState::Compression
        }
    }

    /// [`MemberForce::state_with_tol`] with a default tolerance of
    /// `1e-9` (in the load's force unit).
    #[must_use]
    pub fn state(&self) -> AxialState {
        self.state_with_tol(1e-9)
    }

    /// `true` if the member is in tension (`force > tol`).
    #[must_use]
    pub fn is_tension(&self) -> bool {
        matches!(self.state(), AxialState::Tension)
    }

    /// `true` if the member is in compression (`force < -tol`).
    #[must_use]
    pub fn is_compression(&self) -> bool {
        matches!(self.state(), AxialState::Compression)
    }
}

/// The solved reaction at a supported node, resolved into global
/// `(fx, fy)` components.
///
/// For a [`crate::model::Support::Roller`] the reaction is the single
/// scalar already projected onto its `(fx, fy)` normal, so the same two
/// fields describe both pin and roller reactions uniformly.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Reaction {
    /// Index of the supported node in [`Truss::nodes`].
    pub node: usize,
    /// Horizontal reaction component.
    pub fx: f64,
    /// Vertical reaction component.
    pub fy: f64,
}

/// Internal tag carried from assembly into unpacking, recording whether
/// a reaction column block came from a pin (two columns) or a roller
/// (one column projected onto the stored normal).
#[derive(Debug, Clone, Copy)]
pub(crate) enum ReactionKind {
    /// Pin: two reaction unknowns `(Rx, Ry)` in consecutive columns.
    Pin,
    /// Roller: one reaction unknown acting along the unit normal
    /// `(nx, ny)`.
    Roller {
        /// Normal x-component.
        nx: f64,
        /// Normal y-component.
        ny: f64,
    },
}

/// The complete solution of a determinate truss.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TrussSolution {
    /// One entry per member, in [`Truss::members`] order.
    pub member_forces: Vec<MemberForce>,
    /// One entry per supported node, in node order.
    pub reactions: Vec<Reaction>,
}

impl TrussSolution {
    /// Solved axial force of member index `m`.
    ///
    /// # Panics
    /// Panics if `m` is out of range.
    #[must_use]
    pub fn force(&self, m: usize) -> f64 {
        self.member_forces[m].force
    }

    /// Sum of all reaction `fx` components.
    #[must_use]
    pub fn total_reaction_fx(&self) -> f64 {
        self.reactions.iter().map(|r| r.fx).sum()
    }

    /// Sum of all reaction `fy` components.
    #[must_use]
    pub fn total_reaction_fy(&self) -> f64 {
        self.reactions.iter().map(|r| r.fy).sum()
    }

    /// Global force residual `(ΣFx, ΣFy)` over the whole structure:
    /// every applied joint load plus every support reaction. For a
    /// correctly solved truss both components are zero to within
    /// round-off.
    ///
    /// This is an *independent* equilibrium check — it re-sums the
    /// reactions returned by the solver against the loads read straight
    /// off `truss`, so a non-zero residual would flag an assembly bug.
    #[must_use]
    pub fn global_residual(&self, truss: &Truss) -> (f64, f64) {
        let mut sx = 0.0;
        let mut sy = 0.0;
        for node in &truss.nodes {
            if let Some(load) = node.load {
                sx += load.fx;
                sy += load.fy;
            }
        }
        sx += self.total_reaction_fx();
        sy += self.total_reaction_fy();
        (sx, sy)
    }

    /// Per-joint equilibrium residual at node `i`: the net `(Fx, Fy)`
    /// from all members meeting the joint, its reaction (if any), and
    /// its applied load (if any). Zero to round-off at every joint of a
    /// correctly solved truss — this is the method-of-joints condition
    /// itself, re-evaluated from the solution.
    #[must_use]
    pub fn joint_residual(&self, truss: &Truss, i: usize) -> (f64, f64) {
        let mut fx = 0.0;
        let mut fy = 0.0;

        // Member contributions: a member of tension S meeting joint i
        // pulls the joint toward the far end. Taking the axis from this
        // joint toward the *other* end gives that pull directly, with no
        // sign bookkeeping needed for the a-end vs b-end case.
        for mf in &self.member_forces {
            let other = if mf.a == i {
                mf.b
            } else if mf.b == i {
                mf.a
            } else {
                continue;
            };
            let pi = truss.nodes[i];
            let po = truss.nodes[other];
            let dx = po.x - pi.x;
            let dy = po.y - pi.y;
            let len = (dx * dx + dy * dy).sqrt();
            fx += mf.force * dx / len;
            fy += mf.force * dy / len;
        }

        // Reaction at this joint, if supported.
        for r in &self.reactions {
            if r.node == i {
                fx += r.fx;
                fy += r.fy;
            }
        }

        // Applied load at this joint, if any.
        if let Some(load) = truss.nodes[i].load {
            fx += load.fx;
            fy += load.fy;
        }

        (fx, fy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn axial_state_classification() {
        let mk = |force: f64| MemberForce {
            member: 0,
            a: 0,
            b: 1,
            force,
            length: 1.0,
        };
        assert_eq!(mk(10.0).state(), AxialState::Tension);
        assert_eq!(mk(-10.0).state(), AxialState::Compression);
        assert_eq!(mk(0.0).state(), AxialState::Zero);
        assert!(mk(5.0).is_tension());
        assert!(mk(-5.0).is_compression());
        // Below tolerance counts as zero, not tension.
        assert_eq!(mk(1e-12).state(), AxialState::Zero);
    }

    #[test]
    fn reaction_totals_sum_components() {
        let sol = TrussSolution {
            member_forces: vec![],
            reactions: vec![
                Reaction {
                    node: 0,
                    fx: 1.0,
                    fy: 2.0,
                },
                Reaction {
                    node: 1,
                    fx: -0.5,
                    fy: 3.0,
                },
            ],
        };
        assert!((sol.total_reaction_fx() - 0.5).abs() < 1e-12);
        assert!((sol.total_reaction_fy() - 5.0).abs() < 1e-12);
    }
}
