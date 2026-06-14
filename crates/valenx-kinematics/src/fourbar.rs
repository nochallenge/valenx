//! Planar four-bar linkage: Grashof classification and loop-closure.
//!
//! ## Geometry / sign convention
//!
//! The mechanism is laid out in the plane with the crank pivot `O2` at
//! the origin and the rocker pivot `O4` on the positive x-axis a
//! distance `r1` (the *ground* / frame link) away. Starting from the
//! crank and walking the loop:
//!
//! - `r1` — ground link, from `O2 = (0, 0)` to `O4 = (r1, 0)`.
//! - `r2` — crank (input), from `O2` to moving pin `A`, at angle `θ2`
//!   measured CCW from the +x axis.
//! - `r3` — coupler, from `A` to moving pin `B`.
//! - `r4` — rocker (output / follower), from `O4` to `B`, at angle
//!   `θ4` from the +x axis.
//!
//! The vector loop is `r2 + r3 = r1 + r4`, i.e. the crank pin `A` and
//! the rocker pin `B` are joined by the coupler of fixed length `r3`.
//!
//! ## Solution
//!
//! For a given crank angle `θ2` the crank pin sits at
//! `A = (r2 cosθ2, r2 sinθ2)`. The straight-line distance from `A` to
//! the rocker pivot `O4` (call it the diagonal `e`) is fixed by `θ2`.
//! The rocker pin `B` then lies at the intersection of two circles —
//! one of radius `r3` about `A`, one of radius `r4` about `O4` — which
//! reduces to solving the triangle `A–O4–B` with known sides `e`,
//! `r4`, `r3`. There are (generically) two intersections: the
//! [`Assembly::Open`] and [`Assembly::Crossed`] branches. A real
//! solution exists only when `|r3 − r4| ≤ e ≤ r3 + r4`; otherwise the
//! coupler cannot reach and [`KinematicsError::CannotClose`] is
//! returned.

use serde::{Deserialize, Serialize};

use crate::error::{require_positive, KinematicsError};

/// Grashof mobility class of a planar four-bar, from the link-length
/// inequality `s + l` vs `p + q` where `s` is the shortest link, `l`
/// the longest, and `p`, `q` the remaining two.
///
/// - [`GrashofClass::Crank`] (`s + l < p + q`): at least one link can
///   fully rotate. Which link rotates depends on which is the shortest
///   and which is grounded, but the *class* is Grashof.
/// - [`GrashofClass::Change`] (`s + l == p + q`): the change-point
///   (folding) case; the linkage can pass through a collinear
///   configuration.
/// - [`GrashofClass::DoubleRocker`] (`s + l > p + q`): non-Grashof; no
///   link can make a full revolution relative to any other.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum GrashofClass {
    /// Grashof (`s + l < p + q`): some link fully rotates.
    Crank,
    /// Change-point (`s + l == p + q`): collinear folding case.
    Change,
    /// Non-Grashof (`s + l > p + q`): double-rocker, no full rotation.
    DoubleRocker,
}

/// Which of the two assembly branches (circle–circle intersections)
/// to take when solving the loop. Both satisfy the loop equation; they
/// are mirror images across the crank-pin-to-rocker-pivot diagonal.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Assembly {
    /// The "open" branch (rocker pin on the CCW side of the diagonal).
    Open,
    /// The "crossed" branch (rocker pin on the CW side of the diagonal).
    Crossed,
}

/// A planar four-bar linkage, defined by its four link lengths in a
/// consistent unit. Validated on construction: every length must be a
/// finite value strictly greater than zero.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FourBar {
    /// Ground / frame link length `r1` (`O2`–`O4`).
    pub ground: f64,
    /// Crank (input) link length `r2`.
    pub crank: f64,
    /// Coupler link length `r3`.
    pub coupler: f64,
    /// Rocker (output) link length `r4`.
    pub rocker: f64,
}

/// The solved pose of a four-bar at one crank angle: the moving pin
/// positions and the output-link angles.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FourBarPose {
    /// Crank-pin position `A = O2 + r2·(cosθ2, sinθ2)`.
    pub pin_a: [f64; 2],
    /// Coupler-pin position `B = O4 + r4·(cosθ4, sinθ4)`.
    pub pin_b: [f64; 2],
    /// Coupler absolute angle `θ3` (radians, CCW from +x), the
    /// orientation of the vector from `A` to `B`.
    pub coupler_angle: f64,
    /// Rocker absolute angle `θ4` (radians, CCW from +x), the
    /// orientation of the vector from `O4` to `B`.
    pub rocker_angle: f64,
}

impl FourBar {
    /// Construct and validate a four-bar from its four link lengths.
    ///
    /// # Errors
    /// Returns [`KinematicsError::BadParameter`] if any length is not a
    /// finite, strictly-positive value.
    pub fn new(
        ground: f64,
        crank: f64,
        coupler: f64,
        rocker: f64,
    ) -> Result<Self, KinematicsError> {
        require_positive("ground", ground)?;
        require_positive("crank", crank)?;
        require_positive("coupler", coupler)?;
        require_positive("rocker", rocker)?;
        Ok(Self {
            ground,
            crank,
            coupler,
            rocker,
        })
    }

    /// Classify the linkage under the Grashof criterion.
    ///
    /// Sorts the four lengths, then compares `shortest + longest`
    /// against the sum of the middle two. The comparison uses a small
    /// relative tolerance so that exact change-point linkages (e.g. a
    /// parallelogram, where the sums are equal) are reported as
    /// [`GrashofClass::Change`] rather than flickering between the
    /// strict inequalities through floating-point round-off.
    pub fn grashof_class(&self) -> GrashofClass {
        let mut links = [self.ground, self.crank, self.coupler, self.rocker];
        links.sort_by(|a, b| a.partial_cmp(b).expect("link lengths are finite"));
        let s = links[0];
        let l = links[3];
        let p = links[1];
        let q = links[2];
        let lhs = s + l;
        let rhs = p + q;
        let scale = lhs.abs().max(rhs.abs()).max(1.0);
        let tol = 1e-9 * scale;
        if (lhs - rhs).abs() <= tol {
            GrashofClass::Change
        } else if lhs < rhs {
            GrashofClass::Crank
        } else {
            GrashofClass::DoubleRocker
        }
    }

    /// `true` for a Grashof linkage (some link can fully rotate),
    /// i.e. strictly `s + l < p + q`.
    pub fn is_grashof(&self) -> bool {
        matches!(self.grashof_class(), GrashofClass::Crank)
    }

    /// Solve the loop closure at crank angle `theta2` (radians, CCW
    /// from +x), returning the moving-pin positions and output angles
    /// on the requested [`Assembly`] branch.
    ///
    /// # Errors
    /// Returns [`KinematicsError::CannotClose`] when the coupler cannot
    /// span the gap between the crank pin and the rocker pivot at this
    /// crank angle (the linkage is not assemblable there).
    pub fn solve(&self, theta2: f64, assembly: Assembly) -> Result<FourBarPose, KinematicsError> {
        // Crank pin A.
        let ax = self.crank * theta2.cos();
        let ay = self.crank * theta2.sin();

        // Rocker pivot O4 = (ground, 0). Diagonal from A to O4.
        let dx = self.ground - ax;
        let dy = -ay;
        let e = (dx * dx + dy * dy).sqrt(); // diagonal length

        let r3 = self.coupler;
        let r4 = self.rocker;
        let reach_min = (r3 - r4).abs();
        let reach_max = r3 + r4;

        // Tiny slack so an exactly-grazing configuration (e at the
        // boundary, as at a Grashof toggle position) is treated as
        // closable rather than rejected by round-off.
        let slack = 1e-9 * reach_max.max(1.0);
        if e < reach_min - slack || e > reach_max + slack {
            return Err(KinematicsError::CannotClose {
                crank_rad: theta2,
                diagonal: e,
                reach_min,
                reach_max,
            });
        }

        // Angle of the diagonal (A -> O4) measured at A.
        let diag_angle = dy.atan2(dx);

        // Interior angle at A in the triangle A-O4-B (sides e, r4 and
        // r3) by the law of cosines. Clamp the cosine into [-1, 1] to
        // absorb the boundary round-off allowed by `slack`.
        let cos_a = ((e * e + r3 * r3 - r4 * r4) / (2.0 * e * r3)).clamp(-1.0, 1.0);
        let interior = cos_a.acos();

        // The two branches put B on opposite sides of the diagonal.
        let theta3 = match assembly {
            Assembly::Open => diag_angle - interior,
            Assembly::Crossed => diag_angle + interior,
        };

        // Pin B from the coupler vector off A.
        let bx = ax + r3 * theta3.cos();
        let by = ay + r3 * theta3.sin();

        // Rocker angle is the orientation of O4 -> B.
        let theta4 = (by - 0.0).atan2(bx - self.ground);

        Ok(FourBarPose {
            pin_a: [ax, ay],
            pin_b: [bx, by],
            coupler_angle: theta3,
            rocker_angle: theta4,
        })
    }

    /// Loop-closure residual at a solved pose: the vector
    /// `r2 + r3 − r1 − r4` expressed in the plane. For an exact
    /// solution both components are zero (to floating-point
    /// precision). Useful as a self-check and exercised by the tests.
    pub fn closure_residual(&self, theta2: f64, pose: &FourBarPose) -> [f64; 2] {
        // r2 along the crank.
        let r2 = [self.crank * theta2.cos(), self.crank * theta2.sin()];
        // r3 along the coupler (A -> B).
        let r3 = [pose.pin_b[0] - pose.pin_a[0], pose.pin_b[1] - pose.pin_a[1]];
        // r1 along the ground (O2 -> O4).
        let r1 = [self.ground, 0.0];
        // r4 along the rocker (O4 -> B).
        let r4 = [
            self.rocker * pose.rocker_angle.cos(),
            self.rocker * pose.rocker_angle.sin(),
        ];
        [r2[0] + r3[0] - r1[0] - r4[0], r2[1] + r3[1] - r1[1] - r4[1]]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-9;

    #[test]
    fn rejects_non_positive_links() {
        assert!(FourBar::new(0.0, 1.0, 1.0, 1.0).is_err());
        assert!(FourBar::new(1.0, -1.0, 1.0, 1.0).is_err());
        assert!(FourBar::new(1.0, 1.0, f64::NAN, 1.0).is_err());
        assert!(FourBar::new(1.0, 1.0, 1.0, 1.0).is_ok());
    }

    // --- Grashof classification on known link sets ---

    #[test]
    fn grashof_crank_rocker_known_set() {
        // Classic crank-rocker: ground 4, crank 2, coupler 5, rocker 5.
        // sorted: 2, 4, 5, 5 -> s+l = 7, p+q = 9 -> Grashof.
        let fb = FourBar::new(4.0, 2.0, 5.0, 5.0).unwrap();
        assert_eq!(fb.grashof_class(), GrashofClass::Crank);
        assert!(fb.is_grashof());
    }

    #[test]
    fn non_grashof_double_rocker_known_set() {
        // Lengths 2, 3, 3, 3: s+l = 5, p+q = 6 -> Grashof, so pick a
        // genuinely non-Grashof set: 4, 5, 2, 5 -> sorted 2,4,5,5,
        // that's Grashof too. Use 3, 4, 5, 6 by shortest+longest:
        // sorted 3,4,5,6 -> s+l = 9, p+q = 9 -> change point.
        // A clean double-rocker: 5, 4, 2, 2 -> sorted 2,2,4,5 ->
        // s+l = 7, p+q = 6 -> s+l > p+q -> double rocker.
        let fb = FourBar::new(5.0, 4.0, 2.0, 2.0).unwrap();
        assert_eq!(fb.grashof_class(), GrashofClass::DoubleRocker);
        assert!(!fb.is_grashof());
    }

    #[test]
    fn change_point_parallelogram_and_345_6() {
        // Parallelogram: equal opposite links -> exact change point.
        let par = FourBar::new(4.0, 2.0, 4.0, 2.0).unwrap();
        assert_eq!(par.grashof_class(), GrashofClass::Change);

        // 3,4,5,6: sorted 3,4,5,6 -> s+l = 9 = p+q = 9 -> change.
        let cp = FourBar::new(3.0, 4.0, 5.0, 6.0).unwrap();
        assert_eq!(cp.grashof_class(), GrashofClass::Change);
    }

    #[test]
    fn grashof_independent_of_link_ordering() {
        // Same multiset of lengths in different role assignments must
        // give the same class — the criterion is order-free.
        let a = FourBar::new(4.0, 2.0, 5.0, 5.0).unwrap();
        let b = FourBar::new(5.0, 5.0, 2.0, 4.0).unwrap();
        let c = FourBar::new(2.0, 5.0, 4.0, 5.0).unwrap();
        assert_eq!(a.grashof_class(), b.grashof_class());
        assert_eq!(b.grashof_class(), c.grashof_class());
    }

    // --- Loop closure: vector sum is zero at solved angles ---

    #[test]
    fn loop_closes_for_crank_rocker_sweep_open() {
        let fb = FourBar::new(4.0, 2.0, 5.0, 5.0).unwrap();
        // A crank-rocker's crank fully rotates, so every angle solves.
        let n = 72;
        for i in 0..n {
            let theta2 = 2.0 * PI * (i as f64) / (n as f64);
            let pose = fb
                .solve(theta2, Assembly::Open)
                .unwrap_or_else(|e| panic!("open solve failed at {theta2}: {e}"));
            let res = fb.closure_residual(theta2, &pose);
            assert!(
                res[0].abs() < 1e-9 && res[1].abs() < 1e-9,
                "residual {res:?} too large at theta2 = {theta2}"
            );
        }
    }

    #[test]
    fn loop_closes_for_crank_rocker_sweep_crossed() {
        let fb = FourBar::new(4.0, 2.0, 5.0, 5.0).unwrap();
        let n = 72;
        for i in 0..n {
            let theta2 = 2.0 * PI * (i as f64) / (n as f64);
            let pose = fb.solve(theta2, Assembly::Crossed).unwrap();
            let res = fb.closure_residual(theta2, &pose);
            assert!(
                res[0].abs() < 1e-9 && res[1].abs() < 1e-9,
                "crossed residual {res:?} too large at theta2 = {theta2}"
            );
        }
    }

    #[test]
    fn solved_pin_b_lies_on_both_link_circles() {
        // B must be exactly r3 from A and exactly r4 from O4.
        let fb = FourBar::new(4.0, 2.0, 5.0, 5.0).unwrap();
        let theta2 = 0.7;
        let pose = fb.solve(theta2, Assembly::Open).unwrap();
        let da = ((pose.pin_b[0] - pose.pin_a[0]).powi(2)
            + (pose.pin_b[1] - pose.pin_a[1]).powi(2))
        .sqrt();
        let do4 = ((pose.pin_b[0] - fb.ground).powi(2) + pose.pin_b[1].powi(2)).sqrt();
        assert!(
            (da - fb.coupler).abs() < EPS,
            "|AB| = {da}, want {}",
            fb.coupler
        );
        assert!(
            (do4 - fb.rocker).abs() < EPS,
            "|O4B| = {do4}, want {}",
            fb.rocker
        );
    }

    #[test]
    fn open_and_crossed_branches_differ_then_meet_at_toggle() {
        // Generic interior angle: the two branches put B in different
        // places. At a toggle (interior angle 0) they coincide, but a
        // mid-sweep angle is generically not a toggle.
        let fb = FourBar::new(4.0, 2.0, 5.0, 5.0).unwrap();
        let theta2 = 1.2;
        let open = fb.solve(theta2, Assembly::Open).unwrap();
        let crossed = fb.solve(theta2, Assembly::Crossed).unwrap();
        let gap = ((open.pin_b[0] - crossed.pin_b[0]).powi(2)
            + (open.pin_b[1] - crossed.pin_b[1]).powi(2))
        .sqrt();
        assert!(gap > 1e-3, "branches should differ off-toggle, gap = {gap}");
    }

    #[test]
    fn coupler_and_rocker_angles_are_consistent_with_pins() {
        let fb = FourBar::new(4.0, 2.0, 5.0, 5.0).unwrap();
        let theta2 = 2.3;
        let pose = fb.solve(theta2, Assembly::Open).unwrap();
        // Reconstruct B from the reported angles and compare.
        let bx_from_coupler = pose.pin_a[0] + fb.coupler * pose.coupler_angle.cos();
        let by_from_coupler = pose.pin_a[1] + fb.coupler * pose.coupler_angle.sin();
        assert!((bx_from_coupler - pose.pin_b[0]).abs() < EPS);
        assert!((by_from_coupler - pose.pin_b[1]).abs() < EPS);

        let bx_from_rocker = fb.ground + fb.rocker * pose.rocker_angle.cos();
        let by_from_rocker = fb.rocker * pose.rocker_angle.sin();
        assert!((bx_from_rocker - pose.pin_b[0]).abs() < EPS);
        assert!((by_from_rocker - pose.pin_b[1]).abs() < EPS);
    }

    #[test]
    fn unreachable_configuration_reports_cannot_close() {
        // Non-Grashof double-rocker whose crank cannot reach θ2 = π:
        // ground 5, crank 4, coupler 2, rocker 2. At θ2 = π the crank
        // pin is at (-4, 0); diagonal to O4 = (5,0) is 9, while the
        // coupler+rocker reach is only 4 -> cannot close.
        let fb = FourBar::new(5.0, 4.0, 2.0, 2.0).unwrap();
        let err = fb.solve(PI, Assembly::Open).unwrap_err();
        assert_eq!(err.code(), "kinematics.cannot_close");
        match err {
            KinematicsError::CannotClose {
                diagonal,
                reach_max,
                ..
            } => {
                assert!(
                    diagonal > reach_max,
                    "diag {diagonal} should exceed reach {reach_max}"
                );
            }
            other => panic!("wrong error: {other:?}"),
        }
    }

    #[test]
    fn symmetric_square_linkage_closes_at_right_angle() {
        // All links equal (a rhombus/parallelogram, change-point). At
        // θ2 = π/2 the crank pin is at (0, L). With unit links and
        // ground along x, the open branch should still close cleanly.
        let fb = FourBar::new(1.0, 1.0, 1.0, 1.0).unwrap();
        let pose = fb.solve(PI / 2.0, Assembly::Open).unwrap();
        let res = fb.closure_residual(PI / 2.0, &pose);
        assert!(
            res[0].abs() < 1e-9 && res[1].abs() < 1e-9,
            "residual {res:?}"
        );
    }
}
