//! Method-of-joints assembly and linear solve.
//!
//! The whole truss is solved in one shot as a single linear system
//! `A·x = b`, the matrix form of "ΣF = 0 at every joint". This is the
//! [method of joints] expressed globally instead of joint-by-joint: it
//! solves the identical equations, but lets [`nalgebra`]'s pivoted LU do
//! the elimination and—because it never has to find a joint with only
//! two unknowns to start from—handles any determinate topology, not just
//! the ones that unwind by hand.
//!
//! # The linear system
//!
//! Order the unknowns `x` as `[S₀ … S_{m-1}, R₀ … R_{r-1}]`: first the
//! axial force in each member (tension positive), then the reaction
//! scalars contributed by the supports in node order (two per
//! [`Support::Pin`], one per [`Support::Roller`]).
//!
//! There are two equilibrium equations per node — `ΣFx = 0` and
//! `ΣFy = 0` — so `A` is `(2N) × (m + r)`. For a statically determinate,
//! stable truss `m + r = 2N`, the matrix is square and invertible, and
//!
//! ```text
//!   A · x = b ,   b_node = −(applied load at that node)
//! ```
//!
//! has a unique solution. A member of tension `S` joining node `a` to
//! node `b`, with unit axis `û = (b − a)/‖b − a‖`, *pulls each joint
//! toward the other end*: it adds `+S·û` to the force balance at `a` and
//! `−S·û` at `b`. A pin adds `+Rx`, `+Ry`; a roller adds `+R·n̂` where
//! `n̂` is the unit normal to its slide direction.
//!
//! [method of joints]: https://en.wikipedia.org/wiki/Truss#Analysis

use nalgebra::{DMatrix, DVector};

use crate::error::TrussError;
use crate::model::{Support, Truss, MIN_MEMBER_LENGTH};
use crate::result::{MemberForce, Reaction, ReactionKind, TrussSolution};

/// Assemble and solve the joint-equilibrium system for `truss`,
/// returning the axial member forces and support reactions.
///
/// # Errors
/// - [`TrussError::Empty`] if the truss has no nodes or no members.
/// - [`TrussError::DegenerateMember`] if any member is zero length (only
///   reachable if a [`crate::model::Member`] was constructed by hand and
///   pushed without going through [`Truss::add_member`]).
/// - [`TrussError::ZeroRollerDirection`] if a roller's slide direction is
///   the zero vector.
/// - [`TrussError::NotDeterminate`] if `m + r ≠ 2N` (the count condition
///   for static determinacy fails).
/// - [`TrussError::Singular`] if the square system is rank-deficient
///   (a mechanism or improperly-constrained geometry).
pub fn solve(truss: &Truss) -> Result<TrussSolution, TrussError> {
    let n_nodes = truss.nodes.len();
    if n_nodes == 0 {
        return Err(TrussError::Empty { what: "no nodes" });
    }
    if truss.members.is_empty() {
        return Err(TrussError::Empty { what: "no members" });
    }

    let m = truss.members.len();
    let r = truss.reaction_count();
    let eqs = truss.equation_count(); // 2 * n_nodes
    let unknowns = m + r;

    if unknowns != eqs {
        return Err(TrussError::NotDeterminate {
            unknowns,
            members: m,
            reactions: r,
            equations: eqs,
            nodes: n_nodes,
        });
    }

    let mut a = DMatrix::<f64>::zeros(eqs, unknowns);
    let mut b = DVector::<f64>::zeros(eqs);

    // --- member-force columns ------------------------------------------
    for (k, member) in truss.members.iter().enumerate() {
        let (dx, dy, len) = truss.member_geometry(*member);
        if len < MIN_MEMBER_LENGTH {
            return Err(TrussError::DegenerateMember {
                member: k,
                a: member.a,
                b: member.b,
            });
        }
        let ux = dx / len;
        let uy = dy / len;

        // Tension pulls joint `a` toward `b` (+û) and joint `b`
        // toward `a` (−û).
        a[(2 * member.a, k)] += ux;
        a[(2 * member.a + 1, k)] += uy;
        a[(2 * member.b, k)] -= ux;
        a[(2 * member.b + 1, k)] -= uy;
    }

    // --- reaction columns ----------------------------------------------
    // Reactions are appended after the `m` member columns, in node order.
    let mut col = m;
    // Track which reaction column(s) each supported node owns so the
    // solution can be unpacked back onto nodes afterwards.
    let mut reaction_map: Vec<(usize, ReactionKind, usize)> = Vec::new();
    for (i, node) in truss.nodes.iter().enumerate() {
        let Some(support) = node.support else {
            continue;
        };
        match support {
            Support::Pin => {
                // Rx on the Fx row, Ry on the Fy row.
                a[(2 * i, col)] += 1.0;
                a[(2 * i + 1, col + 1)] += 1.0;
                reaction_map.push((i, ReactionKind::Pin, col));
                col += 2;
            }
            Support::Roller { slide } => {
                let (nx, ny) = roller_normal(i, slide)?;
                a[(2 * i, col)] += nx;
                a[(2 * i + 1, col)] += ny;
                reaction_map.push((i, ReactionKind::Roller { nx, ny }, col));
                col += 1;
            }
        }
    }
    debug_assert_eq!(col, unknowns, "reaction columns must fill exactly");

    // --- right-hand side: −(applied load) ------------------------------
    for (i, node) in truss.nodes.iter().enumerate() {
        if let Some(load) = node.load {
            b[2 * i] -= load.fx;
            b[2 * i + 1] -= load.fy;
        }
    }

    // --- solve ----------------------------------------------------------
    let lu = a.clone().full_piv_lu();
    if !lu.is_invertible() {
        return Err(TrussError::Singular);
    }
    let x = lu.solve(&b).ok_or(TrussError::Singular)?;

    // --- unpack ---------------------------------------------------------
    let member_forces = truss
        .members
        .iter()
        .enumerate()
        .map(|(k, member)| MemberForce {
            member: k,
            a: member.a,
            b: member.b,
            force: x[k],
            length: truss.member_length(k),
        })
        .collect();

    let reactions = reaction_map
        .into_iter()
        .map(|(node, kind, base)| match kind {
            ReactionKind::Pin => Reaction {
                node,
                fx: x[base],
                fy: x[base + 1],
            },
            ReactionKind::Roller { nx, ny } => {
                let scalar = x[base];
                Reaction {
                    node,
                    fx: scalar * nx,
                    fy: scalar * ny,
                }
            }
        })
        .collect();

    Ok(TrussSolution {
        member_forces,
        reactions,
    })
}

/// Unit normal `(nx, ny)` to a roller's slide direction `slide`.
///
/// The normal is the slide vector rotated +90° and normalised, so the
/// roller's single reaction acts perpendicular to the surface it rolls
/// on.
fn roller_normal(node: usize, slide: [f64; 2]) -> Result<(f64, f64), TrussError> {
    let [sx, sy] = slide;
    let mag = (sx * sx + sy * sy).sqrt();
    if mag < MIN_MEMBER_LENGTH {
        return Err(TrussError::ZeroRollerDirection { node });
    }
    // Rotate (sx, sy) by +90°: (-sy, sx); then normalise.
    Ok((-sy / mag, sx / mag))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Load, Member, Node};

    /// Build the canonical loaded triangle used across the tests:
    ///
    /// ```text
    ///            C (2, 3)
    ///           / \
    ///          /   \
    ///   A(0,0)*-----*B(4,0)
    ///   pin        roller(H)
    /// ```
    ///
    /// A vertical down-load `P` is applied at apex C. A is a pin, B is a
    /// horizontal roller (vertical reaction only).
    fn loaded_triangle(p: f64) -> Truss {
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
        t.add_node(
            Node::new(4.0, 0.0)
                .unwrap()
                .with_support(Support::horizontal_roller()),
        );
        t.add_node(
            Node::new(2.0, 3.0)
                .unwrap()
                .with_load(Load::down(p).unwrap()),
        );
        t.add_member(Member::new(0, 1)).unwrap(); // bottom chord AB
        t.add_member(Member::new(1, 2)).unwrap(); // BC
        t.add_member(Member::new(2, 0)).unwrap(); // CA
        t
    }

    #[test]
    fn loaded_triangle_solves_to_analytic_forces() {
        // P = 60: bottom chord +P/3 = +20 (tension), both struts
        // −P·√13/6 = −10·√13 (compression). Members were added in the
        // order AB, BC, CA.
        let sol = solve(&loaded_triangle(60.0)).unwrap();
        let strut = -60.0 * 13.0_f64.sqrt() / 6.0;
        assert!((sol.force(0) - 20.0).abs() < 1e-9);
        assert!((sol.force(1) - strut).abs() < 1e-9);
        assert!((sol.force(2) - strut).abs() < 1e-9);
        assert!(sol.member_forces[0].is_tension());
        assert!(sol.member_forces[1].is_compression());
        // Independent equilibrium check at the loaded apex (node 2).
        let (fx, fy) = sol.joint_residual(&loaded_triangle(60.0), 2);
        assert!(fx.abs() < 1e-9 && fy.abs() < 1e-9);
    }

    #[test]
    fn empty_truss_errors() {
        let t = Truss::new();
        assert!(matches!(solve(&t), Err(TrussError::Empty { .. })));
    }

    #[test]
    fn no_members_errors() {
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
        assert!(matches!(
            solve(&t),
            Err(TrussError::Empty { what: "no members" })
        ));
    }

    #[test]
    fn under_constrained_is_not_determinate() {
        // Triangle but with only a single pin (r = 2). m + r = 3 + 2 = 5,
        // equations = 6 -> not determinate.
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
        t.add_node(Node::new(4.0, 0.0).unwrap());
        t.add_node(Node::new(2.0, 3.0).unwrap());
        t.add_member(Member::new(0, 1)).unwrap();
        t.add_member(Member::new(1, 2)).unwrap();
        t.add_member(Member::new(2, 0)).unwrap();
        match solve(&t) {
            Err(TrussError::NotDeterminate {
                unknowns,
                equations,
                ..
            }) => {
                assert_eq!(unknowns, 5);
                assert_eq!(equations, 6);
            }
            other => panic!("expected NotDeterminate, got {other:?}"),
        }
    }

    #[test]
    fn zero_roller_direction_errors() {
        let mut t = Truss::new();
        t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
        t.add_node(
            Node::new(4.0, 0.0)
                .unwrap()
                .with_support(Support::Roller { slide: [0.0, 0.0] }),
        );
        t.add_node(
            Node::new(2.0, 3.0)
                .unwrap()
                .with_load(Load::down(10.0).unwrap()),
        );
        t.add_member(Member::new(0, 1)).unwrap();
        t.add_member(Member::new(1, 2)).unwrap();
        t.add_member(Member::new(2, 0)).unwrap();
        assert!(matches!(
            solve(&t),
            Err(TrussError::ZeroRollerDirection { node: 1 })
        ));
    }

    #[test]
    fn roller_normal_is_perpendicular_unit() {
        // Horizontal slide (1,0) -> normal (0,1).
        let (nx, ny) = roller_normal(0, [1.0, 0.0]).unwrap();
        assert!((nx - 0.0).abs() < 1e-12);
        assert!((ny - 1.0).abs() < 1e-12);
        // Non-unit slide is normalised.
        let (nx, ny) = roller_normal(0, [0.0, 5.0]).unwrap();
        assert!((nx - (-1.0)).abs() < 1e-12);
        assert!((ny - 0.0).abs() < 1e-12);
    }
}
