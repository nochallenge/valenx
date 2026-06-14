//! Ground-truth integration tests: each truss here is solved by hand
//! first, and the solver's member forces, tension/compression signs,
//! reactions, and global + per-joint equilibrium are checked against
//! those closed-form values.

use valenx_truss::{solve, AxialState, Load, Member, Node, Support, Truss, TrussError};

const EPS: f64 = 1e-9;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < EPS
}

/// A 60-unit-down load at the apex of a symmetric triangle on a pin +
/// horizontal roller. Hand solution:
///   reactions  Ay = By = P/2 = 30, Ax = 0
///   bottom AB  = +P/3 = +20            (tension)
///   struts     = −P·√13 / 6 = −10·√13  (compression, both)
#[test]
fn loaded_triangle_matches_analytic() {
    let p = 60.0;
    let mut t = Truss::new();
    let a = t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
    let b = t.add_node(
        Node::new(4.0, 0.0)
            .unwrap()
            .with_support(Support::horizontal_roller()),
    );
    let c = t.add_node(
        Node::new(2.0, 3.0)
            .unwrap()
            .with_load(Load::down(p).unwrap()),
    );
    let m_ab = t.add_member(Member::new(a, b)).unwrap();
    let m_bc = t.add_member(Member::new(b, c)).unwrap();
    let m_ca = t.add_member(Member::new(c, a)).unwrap();

    let sol = solve(&t).unwrap();

    let strut = -p * 13.0_f64.sqrt() / 6.0; // −10·√13 ≈ −36.0555
    assert!(
        close(sol.force(m_ab), p / 3.0),
        "AB force {}",
        sol.force(m_ab)
    );
    assert!(
        close(sol.force(m_bc), strut),
        "BC force {}",
        sol.force(m_bc)
    );
    assert!(
        close(sol.force(m_ca), strut),
        "CA force {}",
        sol.force(m_ca)
    );

    // Signs: bottom chord tension, both struts compression.
    assert_eq!(sol.member_forces[m_ab].state(), AxialState::Tension);
    assert_eq!(sol.member_forces[m_bc].state(), AxialState::Compression);
    assert_eq!(sol.member_forces[m_ca].state(), AxialState::Compression);

    // Reactions: A pin (0, 30), B roller (0, 30).
    let ra = sol.reactions.iter().find(|r| r.node == a).unwrap();
    let rb = sol.reactions.iter().find(|r| r.node == b).unwrap();
    assert!(close(ra.fx, 0.0) && close(ra.fy, 30.0), "A reaction {ra:?}");
    assert!(close(rb.fx, 0.0) && close(rb.fy, 30.0), "B reaction {rb:?}");

    // Global equilibrium.
    let (rx, ry) = sol.global_residual(&t);
    assert!(
        close(rx, 0.0) && close(ry, 0.0),
        "global residual ({rx}, {ry})"
    );

    // Per-joint equilibrium at every node.
    for i in 0..t.nodes.len() {
        let (fx, fy) = sol.joint_residual(&t, i);
        assert!(
            close(fx, 0.0) && close(fy, 0.0),
            "joint {i} residual ({fx}, {fy})"
        );
    }
}

/// A one-panel cantilever truss anchored to a wall:
///   A(0,0) pin, B(0,3) vertical-roller (horizontal reaction only),
///   C(4,0) carries a 30-unit downward tip load.
/// Members: AC (bottom), BC (diagonal, len 5), AB (vertical, len 3).
/// Hand solution (W = 30):
///   BC = +5W/3 = +50  (tension)
///   AC = −4W/3 = −40  (compression)
///   AB = −W    = −30  (compression)
///   A reaction (4W/3, W) = (40, 30);  B reaction (−4W/3, 0) = (−40, 0)
#[test]
fn cantilever_truss_matches_analytic() {
    let w = 30.0;
    let mut t = Truss::new();
    let a = t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
    let b = t.add_node(
        Node::new(0.0, 3.0)
            .unwrap()
            .with_support(Support::vertical_roller()),
    );
    let c = t.add_node(
        Node::new(4.0, 0.0)
            .unwrap()
            .with_load(Load::down(w).unwrap()),
    );
    let m_ac = t.add_member(Member::new(a, c)).unwrap();
    let m_bc = t.add_member(Member::new(b, c)).unwrap();
    let m_ab = t.add_member(Member::new(a, b)).unwrap();

    let sol = solve(&t).unwrap();

    assert!(
        close(sol.force(m_bc), 5.0 * w / 3.0),
        "BC {}",
        sol.force(m_bc)
    );
    assert!(
        close(sol.force(m_ac), -4.0 * w / 3.0),
        "AC {}",
        sol.force(m_ac)
    );
    assert!(close(sol.force(m_ab), -w), "AB {}", sol.force(m_ab));

    assert!(sol.member_forces[m_bc].is_tension());
    assert!(sol.member_forces[m_ac].is_compression());
    assert!(sol.member_forces[m_ab].is_compression());

    let ra = sol.reactions.iter().find(|r| r.node == a).unwrap();
    let rb = sol.reactions.iter().find(|r| r.node == b).unwrap();
    assert!(
        close(ra.fx, 40.0) && close(ra.fy, 30.0),
        "A reaction {ra:?}"
    );
    assert!(
        close(rb.fx, -40.0) && close(rb.fy, 0.0),
        "B reaction {rb:?}"
    );

    let (rx, ry) = sol.global_residual(&t);
    assert!(
        close(rx, 0.0) && close(ry, 0.0),
        "global residual ({rx}, {ry})"
    );
    for i in 0..t.nodes.len() {
        let (fx, fy) = sol.joint_residual(&t, i);
        assert!(
            close(fx, 0.0) && close(fy, 0.0),
            "joint {i} residual ({fx}, {fy})"
        );
    }
}

/// The vertical member from the apex to the mid-bottom joint of a
/// simply-supported triangle carries no force: at the unloaded
/// mid-bottom joint two collinear bottom members meet one transverse
/// vertical, so the vertical is a textbook **zero-force member**.
///
///         C(2,3)  ← 60 down
///        / | \
///   A(0,0)-D-B(4,0)
///         (2,0)
#[test]
fn zero_force_member_is_detected() {
    let p = 60.0;
    let mut t = Truss::new();
    let a = t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
    let b = t.add_node(
        Node::new(4.0, 0.0)
            .unwrap()
            .with_support(Support::horizontal_roller()),
    );
    let c = t.add_node(
        Node::new(2.0, 3.0)
            .unwrap()
            .with_load(Load::down(p).unwrap()),
    );
    let d = t.add_node(Node::new(2.0, 0.0).unwrap());

    let m_ad = t.add_member(Member::new(a, d)).unwrap();
    let m_db = t.add_member(Member::new(d, b)).unwrap();
    let m_ac = t.add_member(Member::new(a, c)).unwrap();
    let m_cb = t.add_member(Member::new(c, b)).unwrap();
    let m_cd = t.add_member(Member::new(c, d)).unwrap(); // the vertical

    assert!(t.is_count_determinate());
    let sol = solve(&t).unwrap();

    // The vertical CD is the zero-force member.
    assert_eq!(sol.member_forces[m_cd].state(), AxialState::Zero);
    assert!(close(sol.force(m_cd), 0.0), "CD force {}", sol.force(m_cd));

    // With CD = 0 the structure is the plain triangle: bottom chords
    // each carry +P/3, struts −10·√13.
    let strut = -p * 13.0_f64.sqrt() / 6.0;
    assert!(close(sol.force(m_ad), p / 3.0));
    assert!(close(sol.force(m_db), p / 3.0));
    assert!(close(sol.force(m_ac), strut));
    assert!(close(sol.force(m_cb), strut));

    let (rx, ry) = sol.global_residual(&t);
    assert!(close(rx, 0.0) && close(ry, 0.0));
}

/// Doubling the load doubles every member force and reaction: the system
/// is linear, so the solution scales exactly with the right-hand side.
#[test]
fn solution_scales_linearly_with_load() {
    fn triangle(p: f64) -> Truss {
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
        t.add_member(Member::new(0, 1)).unwrap();
        t.add_member(Member::new(1, 2)).unwrap();
        t.add_member(Member::new(2, 0)).unwrap();
        t
    }

    let s1 = solve(&triangle(10.0)).unwrap();
    let s2 = solve(&triangle(20.0)).unwrap();
    for (a, b) in s1.member_forces.iter().zip(&s2.member_forces) {
        assert!(close(2.0 * a.force, b.force), "{} vs {}", a.force, b.force);
    }
}

/// A horizontal load on the loaded-triangle apex breaks the left/right
/// symmetry; the solver must still satisfy global equilibrium exactly.
#[test]
fn asymmetric_load_still_balances() {
    let mut t = Truss::new();
    let a = t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
    let b = t.add_node(
        Node::new(4.0, 0.0)
            .unwrap()
            .with_support(Support::horizontal_roller()),
    );
    // Combined horizontal + vertical apex load.
    let c = t.add_node(
        Node::new(2.0, 3.0)
            .unwrap()
            .with_load(Load::new(25.0, -60.0).unwrap()),
    );
    t.add_member(Member::new(a, b)).unwrap();
    t.add_member(Member::new(b, c)).unwrap();
    t.add_member(Member::new(c, a)).unwrap();

    let sol = solve(&t).unwrap();

    // Only the pin can take horizontal load (the roller reacts purely
    // vertically), so ΣFx = 0 forces the pin's fx to cancel the applied
    // +25 exactly.
    let ra = sol.reactions.iter().find(|r| r.node == a).unwrap();
    assert!(close(ra.fx, -25.0), "pin fx {}", ra.fx);

    let (rx, ry) = sol.global_residual(&t);
    assert!(close(rx, 0.0) && close(ry, 0.0), "residual ({rx}, {ry})");
    for i in 0..t.nodes.len() {
        let (fx, fy) = sol.joint_residual(&t, i);
        assert!(close(fx, 0.0) && close(fy, 0.0), "joint {i} ({fx}, {fy})");
    }
}

/// A square with a single diagonal, supported by a pin and a roller and
/// loaded at a free top corner, is statically determinate
/// (m = 5, r = 3, N = 4 ⇒ m + r = 2N = 8) and solves cleanly with global
/// equilibrium satisfied.
#[test]
fn square_with_diagonal_is_determinate_and_balances() {
    // Corners of a unit-ish square:
    //   0(0,0) pin     1(4,0) roller(H)
    //   3(0,4)         2(4,4) ← load
    let mut t = Truss::new();
    let n0 = t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
    let n1 = t.add_node(
        Node::new(4.0, 0.0)
            .unwrap()
            .with_support(Support::horizontal_roller()),
    );
    let n2 = t.add_node(
        Node::new(4.0, 4.0)
            .unwrap()
            .with_load(Load::down(40.0).unwrap()),
    );
    let n3 = t.add_node(Node::new(0.0, 4.0).unwrap());

    t.add_member(Member::new(n0, n1)).unwrap(); // bottom
    t.add_member(Member::new(n1, n2)).unwrap(); // right
    t.add_member(Member::new(n2, n3)).unwrap(); // top
    t.add_member(Member::new(n3, n0)).unwrap(); // left
    t.add_member(Member::new(n0, n2)).unwrap(); // diagonal

    assert!(t.is_count_determinate());
    let sol = solve(&t).unwrap();

    let (rx, ry) = sol.global_residual(&t);
    assert!(close(rx, 0.0) && close(ry, 0.0), "residual ({rx}, {ry})");
    for i in 0..t.nodes.len() {
        let (fx, fy) = sol.joint_residual(&t, i);
        assert!(close(fx, 0.0) && close(fy, 0.0), "joint {i} ({fx}, {fy})");
    }
    // Total vertical reaction must carry the whole 40-unit load.
    assert!(
        close(sol.total_reaction_fy(), 40.0),
        "ΣRy {}",
        sol.total_reaction_fy()
    );
}

/// A count-balanced but geometrically unstable truss (a *mechanism*) is
/// rejected as singular rather than returning a bogus force state. Here
/// all three nodes are collinear, so the structure can fold.
#[test]
fn collinear_mechanism_is_singular() {
    let mut t = Truss::new();
    t.add_node(Node::new(0.0, 0.0).unwrap().with_support(Support::Pin));
    t.add_node(
        Node::new(4.0, 0.0)
            .unwrap()
            .with_support(Support::horizontal_roller()),
    );
    // Third node ON the line AB -> degenerate (collinear) triangle.
    t.add_node(
        Node::new(2.0, 0.0)
            .unwrap()
            .with_load(Load::down(10.0).unwrap()),
    );
    t.add_member(Member::new(0, 1)).unwrap();
    t.add_member(Member::new(1, 2)).unwrap();
    t.add_member(Member::new(2, 0)).unwrap();

    assert!(t.is_count_determinate()); // counting passes …
    assert!(matches!(solve(&t), Err(TrussError::Singular))); // … but rank fails.
}
