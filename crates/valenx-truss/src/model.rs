//! Geometry and topology of a planar pin-jointed truss.
//!
//! A [`Truss`] is an ordered list of [`Node`]s (each optionally carrying
//! a [`Support`] and a point [`Load`]) plus a list of two-force
//! [`Member`]s joining node pairs. Everything here is pure data with
//! *validated* mutators: a coordinate is rejected if non-finite, a
//! member is rejected if its endpoints are out of range or coincident,
//! and a roller's slide direction is rejected if it is the zero vector.
//! Solving is done separately in [`crate::solver`].
//!
//! # Sign and frame conventions
//!
//! - The plane is the right-handed `(x, y)` plane; `+x` points right,
//!   `+y` points up. Gravity, if you model it, is a downward (`-y`)
//!   load you apply yourself — members are weightless.
//! - A member axial force is **tension-positive**: a positive solved
//!   force means the member is stretched and pulls its two end joints
//!   *toward* each other; a negative force means compression.
//! - A [`Support::Pin`] removes both translational degrees of freedom
//!   (two reaction unknowns, `Rx` and `Ry`). A [`Support::Roller`]
//!   removes one: it is free to slide along its `slide` direction and
//!   reacts only along the perpendicular *normal*, contributing a single
//!   reaction unknown.

use serde::{Deserialize, Serialize};

use crate::error::TrussError;

/// Members shorter than this are treated as degenerate (zero length, no
/// defined axis). A tiny absolute geometric tolerance, well below any
/// physically meaningful member length in the unit systems the crate
/// targets, yet large enough to reject true coincidences robustly.
pub(crate) const MIN_MEMBER_LENGTH: f64 = 1e-12;

/// A pinned joint of the truss at a fixed `(x, y)` location.
///
/// Each node may carry an optional [`Support`] (otherwise it is free)
/// and an optional point [`Load`] (otherwise unloaded).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Node {
    /// Horizontal coordinate.
    pub x: f64,
    /// Vertical coordinate.
    pub y: f64,
    /// Boundary support at this joint, if any.
    pub support: Option<Support>,
    /// Externally applied point load at this joint, if any.
    pub load: Option<Load>,
}

impl Node {
    /// A free (unsupported, unloaded) node at `(x, y)`.
    ///
    /// # Errors
    /// [`TrussError::NonFinite`] if either coordinate is not finite.
    pub fn new(x: f64, y: f64) -> Result<Self, TrussError> {
        finite("node.x", x)?;
        finite("node.y", y)?;
        Ok(Self {
            x,
            y,
            support: None,
            load: None,
        })
    }

    /// Return a copy of this node with the given [`Support`] attached.
    #[must_use]
    pub fn with_support(mut self, support: Support) -> Self {
        self.support = Some(support);
        self
    }

    /// Return a copy of this node with the given [`Load`] attached.
    #[must_use]
    pub fn with_load(mut self, load: Load) -> Self {
        self.load = Some(load);
        self
    }
}

/// A boundary support fixing a joint against motion.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Support {
    /// A pin (hinge): both `x` and `y` translations are fixed. Supplies
    /// two reaction components `(Rx, Ry)` and no moment.
    Pin,
    /// A roller / link: the joint is free to slide along the unit
    /// direction `slide` and reacts only along the perpendicular
    /// normal. Supplies a single scalar reaction.
    Roller {
        /// Direction (need not be unit length; it is normalised) the
        /// joint is free to translate along.
        slide: [f64; 2],
    },
}

impl Support {
    /// A roller free to slide horizontally (reacts vertically) — the
    /// textbook "roller on a horizontal floor".
    #[must_use]
    pub fn horizontal_roller() -> Self {
        Support::Roller { slide: [1.0, 0.0] }
    }

    /// A roller free to slide vertically (reacts horizontally) — a
    /// roller bearing against a vertical wall.
    #[must_use]
    pub fn vertical_roller() -> Self {
        Support::Roller { slide: [0.0, 1.0] }
    }

    /// Number of scalar reaction unknowns this support introduces:
    /// `2` for a [`Support::Pin`], `1` for a [`Support::Roller`].
    #[must_use]
    pub fn reaction_count(&self) -> usize {
        match self {
            Support::Pin => 2,
            Support::Roller { .. } => 1,
        }
    }
}

/// An externally applied point load at a joint, in global components.
///
/// Components are *forces*, in whatever consistent unit system the
/// caller uses (N, kN, lbf, …). `fy` is positive *upward*.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Load {
    /// Horizontal force component (positive `+x`).
    pub fx: f64,
    /// Vertical force component (positive `+y`).
    pub fy: f64,
}

impl Load {
    /// A load with the given `(fx, fy)` components.
    ///
    /// # Errors
    /// [`TrussError::NonFinite`] if either component is not finite.
    pub fn new(fx: f64, fy: f64) -> Result<Self, TrussError> {
        finite("load.fx", fx)?;
        finite("load.fy", fy)?;
        Ok(Self { fx, fy })
    }

    /// A purely downward load of magnitude `mag` (i.e. `fy = -mag`).
    ///
    /// # Errors
    /// [`TrussError::NonFinite`] if `mag` is not finite.
    pub fn down(mag: f64) -> Result<Self, TrussError> {
        finite("load.mag", mag)?;
        Ok(Self { fx: 0.0, fy: -mag })
    }
}

/// A straight, weightless, two-force member joining two nodes by their
/// indices into [`Truss::nodes`].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Member {
    /// Index of the first end node.
    pub a: usize,
    /// Index of the second end node.
    pub b: usize,
}

impl Member {
    /// A member between node indices `a` and `b`.
    ///
    /// Range and degeneracy are not checked here (the node coordinates
    /// are not yet known); they are validated by
    /// [`Truss::add_member`].
    #[must_use]
    pub fn new(a: usize, b: usize) -> Self {
        Self { a, b }
    }
}

/// A complete planar pin-jointed truss: nodes, members, and (carried on
/// the nodes) supports and loads.
///
/// Build it incrementally — [`Truss::add_node`] returns the new node's
/// index for use in [`Truss::add_member`] — then hand it to
/// [`crate::solver::solve`].
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Truss {
    /// The joints, in insertion order; an index into this vector is a
    /// node id.
    pub nodes: Vec<Node>,
    /// The members joining node pairs.
    pub members: Vec<Member>,
}

impl Truss {
    /// An empty truss.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a node and return its index (its id for member wiring).
    pub fn add_node(&mut self, node: Node) -> usize {
        let id = self.nodes.len();
        self.nodes.push(node);
        id
    }

    /// Append a member after validating its endpoints.
    ///
    /// # Errors
    /// - [`TrussError::NodeIndexOutOfRange`] if `a` or `b` is not a
    ///   valid node index.
    /// - [`TrussError::DegenerateMember`] if `a == b` or the two nodes
    ///   are coincident (zero length, so the axis is undefined).
    pub fn add_member(&mut self, member: Member) -> Result<usize, TrussError> {
        let n = self.nodes.len();
        if member.a >= n {
            return Err(TrussError::NodeIndexOutOfRange {
                index: member.a,
                count: n,
            });
        }
        if member.b >= n {
            return Err(TrussError::NodeIndexOutOfRange {
                index: member.b,
                count: n,
            });
        }
        let idx = self.members.len();
        let (_, _, len) = self.member_geometry(member);
        // Coincident endpoints (a == b is the exact-zero sub-case) leave
        // the member axis undefined. Coordinates are validated finite, so
        // `len` is finite and non-negative; reject anything that rounds to
        // zero against a tight geometric tolerance.
        if member.a == member.b || len < MIN_MEMBER_LENGTH {
            return Err(TrussError::DegenerateMember {
                member: idx,
                a: member.a,
                b: member.b,
            });
        }
        self.members.push(member);
        Ok(idx)
    }

    /// Total number of scalar reaction unknowns across all supported
    /// nodes (`2` per pin, `1` per roller).
    #[must_use]
    pub fn reaction_count(&self) -> usize {
        self.nodes
            .iter()
            .filter_map(|node| node.support.as_ref())
            .map(Support::reaction_count)
            .sum()
    }

    /// Number of joint-equilibrium equations, `2 · N` (two per node).
    #[must_use]
    pub fn equation_count(&self) -> usize {
        2 * self.nodes.len()
    }

    /// Number of unknowns `m + r`: one axial force per member plus the
    /// total reaction components.
    #[must_use]
    pub fn unknown_count(&self) -> usize {
        self.members.len() + self.reaction_count()
    }

    /// Whether the unknown count equals the equation count — the
    /// necessary count condition for static determinacy (`m + r = 2·N`).
    ///
    /// This is the *counting* check only; a count-balanced truss can
    /// still be a geometric mechanism, which the solver catches as
    /// [`TrussError::Singular`].
    #[must_use]
    pub fn is_count_determinate(&self) -> bool {
        self.unknown_count() == self.equation_count()
    }

    /// Geometry of a member: `(dx, dy, length)` from node `a` to node
    /// `b`. The direction cosines are `(dx/length, dy/length)`.
    pub(crate) fn member_geometry(&self, m: Member) -> (f64, f64, f64) {
        let pa = self.nodes[m.a];
        let pb = self.nodes[m.b];
        let dx = pb.x - pa.x;
        let dy = pb.y - pa.y;
        (dx, dy, (dx * dx + dy * dy).sqrt())
    }

    /// Euclidean length of member index `m`.
    ///
    /// # Panics
    /// Panics if `m` is out of range. (Members added via
    /// [`Truss::add_member`] are always in range and non-degenerate.)
    #[must_use]
    pub fn member_length(&self, m: usize) -> f64 {
        let (_, _, len) = self.member_geometry(self.members[m]);
        len
    }
}

/// Reject a non-finite scalar with a named [`TrussError::NonFinite`].
pub(crate) fn finite(what: &'static str, value: f64) -> Result<(), TrussError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(TrussError::NonFinite { what, value })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn node_rejects_nan_coordinate() {
        let err = Node::new(f64::NAN, 0.0).unwrap_err();
        assert_eq!(err.code(), "truss.non_finite");
        assert_eq!(err.category(), crate::error::ErrorCategory::Input);
    }

    #[test]
    fn node_rejects_infinite_coordinate() {
        assert!(matches!(
            Node::new(0.0, f64::INFINITY),
            Err(TrussError::NonFinite { what: "node.y", .. })
        ));
    }

    #[test]
    fn load_down_points_negative_y() {
        let l = Load::down(50.0).unwrap();
        assert!((l.fx - 0.0).abs() < 1e-12);
        assert!((l.fy - (-50.0)).abs() < 1e-12);
    }

    #[test]
    fn support_reaction_counts() {
        assert_eq!(Support::Pin.reaction_count(), 2);
        assert_eq!(Support::horizontal_roller().reaction_count(), 1);
        assert_eq!(Support::vertical_roller().reaction_count(), 1);
    }

    #[test]
    fn member_out_of_range_rejected() {
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap());
        let err = t.add_member(Member::new(0, 5)).unwrap_err();
        assert!(matches!(
            err,
            TrussError::NodeIndexOutOfRange { index: 5, count: 1 }
        ));
    }

    #[test]
    fn member_self_loop_rejected() {
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap());
        assert!(matches!(
            t.add_member(Member::new(0, 0)),
            Err(TrussError::DegenerateMember { a: 0, b: 0, .. })
        ));
    }

    #[test]
    fn coincident_nodes_make_degenerate_member() {
        let mut t = Truss::new();
        t.add_node(Node::new(1.0, 1.0).unwrap());
        t.add_node(Node::new(1.0, 1.0).unwrap());
        assert!(matches!(
            t.add_member(Member::new(0, 1)),
            Err(TrussError::DegenerateMember { .. })
        ));
    }

    #[test]
    fn member_length_is_euclidean() {
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap());
        t.add_node(Node::new(3.0, 4.0).unwrap());
        t.add_member(Member::new(0, 1)).unwrap();
        assert!((t.member_length(0) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn counting_determinacy() {
        // Single triangle: 3 nodes, 3 members. Pin + horizontal roller
        // => r = 3. m + r = 6 = 2*3 nodes. Determinate.
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
        t.add_node(
            Node::new(4.0, 0.0)
                .unwrap()
                .with_support(Support::horizontal_roller()),
        );
        t.add_node(Node::new(2.0, 3.0).unwrap());
        t.add_member(Member::new(0, 1)).unwrap();
        t.add_member(Member::new(1, 2)).unwrap();
        t.add_member(Member::new(2, 0)).unwrap();
        assert_eq!(t.unknown_count(), 6);
        assert_eq!(t.equation_count(), 6);
        assert!(t.is_count_determinate());
    }
}
