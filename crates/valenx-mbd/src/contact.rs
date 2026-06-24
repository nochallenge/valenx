//! Penalty (compliant) **contact dynamics** with regularized **Coulomb
//! friction** for the planar solver.
//!
//! This is a third, complementary model alongside the constrained-DAE
//! [`crate::System`] and the articulated-body [`crate::aba`]: instead of an
//! *exact* non-penetration constraint (a Lagrange multiplier / complementarity
//! problem), a contact is treated as a stiff **spring–damper** that switches on
//! only while a body point has penetrated a rigid half-space (the ground /
//! plane). This is the classic *penalty* (regularized) contact formulation —
//! the same family a game-physics or preliminary-design code uses — and it
//! composes directly with the planar engine, which already assembles applied
//! forces as generalised wrenches `(Fx, Fy, τ)` per body.
//!
//! ## The model
//!
//! A [`Plane`] is an oriented line in the plane: an `outside` half-space (where
//! the inward unit normal `n̂` points) is collision-free, and the opposite side
//! is "inside the ground". For a point at `p` with velocity `ṗ`:
//!
//! - **Penetration depth** `d = n̂·(p₀ − p)` where `p₀` is any point on the
//!   plane; `d > 0` means the point has crossed to the inside by `d` metres.
//! - **Normal force** (only while `d > 0`, and only *pushing out* — never
//!   pulling):
//!   `fₙ = max(0, kₙ·d + cₙ·ḋ)`, with `ḋ = −n̂·ṗ` the closing speed
//!   (positive while sinking deeper). `kₙ` is the contact stiffness (N/m) and
//!   `cₙ` the contact damping (N·s/m). The force acts along `+n̂`. Clamping at
//!   zero keeps the contact *unilateral* — the damper can never yank a
//!   separating point back.
//! - **Tangential / friction force**: a regularized Coulomb law bounded by the
//!   friction cone `|f_t| ≤ μ·fₙ`. Let `v_t` be the tangential component of
//!   `ṗ`. The ideal Coulomb force is `−μ·fₙ·sign(v_t)`; the discontinuity at
//!   `v_t = 0` causes stick–slip chatter, so it is smoothed to
//!   `f_t = −μ·fₙ·tanh(|v_t| / v_eps)·t̂`. For `|v_t| ≫ v_eps` this is the
//!   full kinetic friction `μ·fₙ` opposing motion; for `|v_t| ≲ v_eps` it
//!   ramps linearly through zero (a stiff viscous "stick" region) instead of
//!   flipping sign every step.
//!
//! The result, [`contact_force`], is the **total** contact force `f = fₙ·n̂ +
//! f_t·t̂` on the point. [`body_contact_wrench`] then sums such forces over a
//! body's contact points and reduces them to a generalised wrench
//! `(force, torque about the centre of mass)` ready to drop into the planar
//! solver's applied-force assembly.
//!
//! ## Honest scope
//!
//! Penalty contact is *soft*: there is always a small steady penetration (a
//! resting body sinks `m·g/kₙ` — pinned exactly by the analytic tests), and a
//! very stiff `kₙ` needs a small time step for the explicit integrator to stay
//! stable.
//! This is the documented trade-off of the method (no constraint drift bookkeeping,
//! trivially composable, but compliant and step-size sensitive) — it is a real
//! v1 of contact for the planar engine, not a substitute for a hard
//! complementarity / LCP contact solver. Plane-vs-point only (the half-space is
//! the canonical ground); body-vs-body contact geometry is a documented next
//! step.

use nalgebra::Vector2;

/// 2-D perpendicular (90° CCW): `perp((x, y)) = (−y, x)`. Kept local to this
/// module (the crate root has an identical private helper).
fn perp(v: Vector2<f64>) -> Vector2<f64> {
    Vector2::new(-v.y, v.x)
}

/// 2-D scalar cross product `a × b = aₓ·b_y − a_y·bₓ` (the out-of-plane
/// component). The moment of a planar force `f` applied at arm `r` about the
/// origin is `r × f`.
fn cross2(a: Vector2<f64>, b: Vector2<f64>) -> f64 {
    a.x * b.y - a.y * b.x
}

/// An oriented contact plane (a line in 2-D bounding a rigid half-space).
///
/// `normal` is the **inward** unit normal — it points *out of* the ground into
/// the free half-space, i.e. toward the side a body should stay on. A point is
/// *penetrating* when it lies on the far side of the plane from `normal`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    /// A point lying on the plane (m).
    pub point: Vector2<f64>,
    /// Inward unit normal (points into the collision-free half-space). Need not
    /// be supplied unit — [`Plane::new`] normalises it.
    pub normal: Vector2<f64>,
}

impl Plane {
    /// A plane through `point` with inward `normal` (normalised here). The
    /// normal must be non-zero and finite.
    ///
    /// # Panics
    /// Panics if `normal` is the zero vector or contains a non-finite
    /// component — a degenerate plane has no defined penetration direction.
    pub fn new(point: Vector2<f64>, normal: Vector2<f64>) -> Self {
        assert!(
            point.x.is_finite() && point.y.is_finite(),
            "Plane::new: point must be finite, got {point:?}"
        );
        let len = normal.norm();
        assert!(
            len.is_finite() && len > 0.0,
            "Plane::new: normal must be non-zero and finite, got {normal:?}"
        );
        Self {
            point,
            normal: normal / len,
        }
    }

    /// The flat **ground** plane `y = height`, with the free half-space above
    /// it (inward normal `+y`). The common case.
    pub fn ground(height: f64) -> Self {
        Self::new(Vector2::new(0.0, height), Vector2::new(0.0, 1.0))
    }

    /// Signed **penetration depth** of `p`: `n̂·(p₀ − p)`. Positive when `p`
    /// has crossed into the ground (the penetrating side); zero on the surface;
    /// negative when `p` is safely in the free half-space.
    pub fn penetration(&self, p: Vector2<f64>) -> f64 {
        self.normal.dot(&(self.point - p))
    }
}

/// Penalty-contact parameters: normal spring/damper gains, Coulomb friction
/// coefficient, and the friction-regularization velocity scale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactParams {
    /// Normal contact stiffness `kₙ` (N/m). Must be `> 0`.
    pub normal_stiffness: f64,
    /// Normal contact damping `cₙ` (N·s/m). Must be `≥ 0` (≥ 0, not > 0, so a
    /// purely elastic contact is allowed).
    pub normal_damping: f64,
    /// Coulomb friction coefficient `μ` (dimensionless). Must be `≥ 0`
    /// (`0` = frictionless).
    pub friction: f64,
    /// Tangential-velocity regularization scale `v_eps` (m/s): below this the
    /// friction force ramps smoothly through zero instead of switching sign
    /// (kills stick–slip chatter). Must be `> 0`.
    pub friction_vel_scale: f64,
}

impl ContactParams {
    /// Construct and **validate** a parameter set.
    ///
    /// # Panics
    /// Panics (fail-loud) on a physically invalid configuration:
    /// `normal_stiffness ≤ 0`, `normal_damping < 0`, `friction < 0`,
    /// `friction_vel_scale ≤ 0`, or any non-finite (NaN/∞) field.
    pub fn new(
        normal_stiffness: f64,
        normal_damping: f64,
        friction: f64,
        friction_vel_scale: f64,
    ) -> Self {
        assert!(
            normal_stiffness.is_finite() && normal_stiffness > 0.0,
            "ContactParams: normal_stiffness must be finite and > 0, got {normal_stiffness}"
        );
        assert!(
            normal_damping.is_finite() && normal_damping >= 0.0,
            "ContactParams: normal_damping must be finite and >= 0, got {normal_damping}"
        );
        assert!(
            friction.is_finite() && friction >= 0.0,
            "ContactParams: friction (mu) must be finite and >= 0, got {friction}"
        );
        assert!(
            friction_vel_scale.is_finite() && friction_vel_scale > 0.0,
            "ContactParams: friction_vel_scale must be finite and > 0, got {friction_vel_scale}"
        );
        Self {
            normal_stiffness,
            normal_damping,
            friction,
            friction_vel_scale,
        }
    }
}

/// Total penalty-contact **force** on a point at `position` moving with
/// `velocity`, against `plane`, under `params`.
///
/// Returns `fₙ·n̂ + f_t·t̂`:
/// - `fₙ = max(0, kₙ·d + cₙ·ḋ)` while penetration `d > 0` (else the point is
///   clear of the ground and the force is exactly zero), with closing speed
///   `ḋ = −n̂·velocity`. The `max(0, …)` keeps it compressive — a contact can
///   push out but never pull a point back.
/// - `f_t = −μ·fₙ·tanh(|v_t| / v_eps)`, opposing the tangential velocity
///   `v_t`; bounded by the friction cone `|f_t| ≤ μ·fₙ` by construction
///   (`|tanh| < 1`), and regularized at `v_t ≈ 0` so it does not chatter.
///
/// The friction direction is only formed when `|v_t|` is meaningfully non-zero
/// (guards the divide); at `v_t ≈ 0` the friction force is zero.
pub fn contact_force(
    position: Vector2<f64>,
    velocity: Vector2<f64>,
    plane: &Plane,
    params: &ContactParams,
) -> Vector2<f64> {
    let n = plane.normal; // already unit (Plane::new normalises)
    let depth = plane.penetration(position);
    // Not penetrating (depth ≤ 0), or a non-finite depth from a bad/NaN input
    // position → no contact force at all. Written `is_finite() || depth <= 0`
    // (not `!(depth > 0.0)`) so NaN is rejected explicitly and clippy is happy.
    if !depth.is_finite() || depth <= 0.0 {
        return Vector2::zeros();
    }

    // Closing speed: rate at which penetration is increasing.
    // ḋ = d/dt [ n̂·(p₀ − p) ] = −n̂·ṗ.
    let closing_speed = -n.dot(&velocity);

    // Normal penalty force, clamped to be compressive (no pull).
    let f_n = (params.normal_stiffness * depth + params.normal_damping * closing_speed).max(0.0);
    if f_n == 0.0 {
        // Penetrating but the damper cancels the spring (a fast-separating
        // point at shallow depth): no net contact force this instant.
        return Vector2::zeros();
    }

    let normal_force = n * f_n;

    if params.friction == 0.0 {
        return normal_force;
    }

    // Tangential velocity = velocity with its normal component removed.
    let v_normal = n.dot(&velocity);
    let v_tang_vec = velocity - n * v_normal;
    let v_tang_mag = v_tang_vec.norm();

    // Guard the divide: only form a tangent direction when there is a
    // meaningful sliding speed. `tanh` regularizes the magnitude so that for
    // |v_t| ≲ v_eps the force ramps linearly through zero (stick) and for
    // |v_t| ≫ v_eps it saturates at kinetic friction μ·fₙ (slip).
    let friction_force = if v_tang_mag > 1e-12 {
        let t_hat = v_tang_vec / v_tang_mag;
        let mag = params.friction * f_n * (v_tang_mag / params.friction_vel_scale).tanh();
        -t_hat * mag
    } else {
        Vector2::zeros()
    };

    normal_force + friction_force
}

/// A contact point fixed in a body's **local** frame, paired with the plane it
/// may collide against. Feed a slice of these to [`body_contact_wrench`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyContact {
    /// Contact point in the body's local frame (m), e.g. a corner of a block.
    pub local: Vector2<f64>,
    /// The plane this point collides with.
    pub plane: Plane,
    /// Penalty parameters for this contact.
    pub params: ContactParams,
}

/// The minimal rigid-body state [`body_contact_wrench`] needs — a structural
/// subset of [`crate::Body`]. Centre-of-mass position/orientation and their
/// rates are enough to place every body-fixed contact point and find its
/// velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactBodyState {
    /// Centre-of-mass position (m).
    pub pos: Vector2<f64>,
    /// Orientation (rad).
    pub angle: f64,
    /// Centre-of-mass velocity (m/s).
    pub vel: Vector2<f64>,
    /// Angular velocity (rad/s).
    pub omega: f64,
}

impl From<&crate::Body> for ContactBodyState {
    fn from(b: &crate::Body) -> Self {
        Self {
            pos: b.pos,
            angle: b.angle,
            vel: b.vel,
            omega: b.omega,
        }
    }
}

/// Generalised contact wrench on a body: a planar force and the torque it
/// produces about the body's centre of mass.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Wrench {
    /// Net contact force on the body (N), world frame.
    pub force: Vector2<f64>,
    /// Net contact torque about the centre of mass (N·m); the out-of-plane
    /// component, sign by the right-hand rule (CCW positive).
    pub torque: f64,
}

impl Wrench {
    /// The zero wrench.
    pub fn zero() -> Self {
        Self {
            force: Vector2::zeros(),
            torque: 0.0,
        }
    }
}

/// Sum the penalty-contact forces over a body's contact points and reduce them
/// to a single [`Wrench`] `(force, torque about CoM)` on the body.
///
/// For each contact the body-local point is carried to the world by the body's
/// pose `R(θ)`, its world velocity is `v_cm + ω × r` (planar: `v_cm + ω·perp(r)`),
/// [`contact_force`] gives the force there, and the force is accumulated along
/// with its moment `r × f` about the centre of mass. The returned force/torque
/// drop straight into the planar solver's generalised force `(Fx, Fy, τ)`.
pub fn body_contact_wrench(body: &ContactBodyState, contacts: &[BodyContact]) -> Wrench {
    let (s, c) = body.angle.sin_cos();
    let mut total = Wrench::zero();
    for contact in contacts {
        // World arm from the CoM to the contact point: r = R(θ)·local.
        let l = contact.local;
        let arm = Vector2::new(c * l.x - s * l.y, s * l.x + c * l.y);
        let point_world = body.pos + arm;
        // Velocity of that material point: v_cm + ω × r (planar cross product).
        let point_vel = body.vel + body.omega * perp(arm);
        let f = contact_force(point_world, point_vel, &contact.plane, &contact.params);
        total.force += f;
        total.torque += cross2(arm, f);
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Body;

    const G: f64 = 9.81;

    /// (1) A point mass **resting** on the ground reaches equilibrium with the
    /// normal force equal to its weight and a steady penetration `m·g/kₙ` —
    /// matched to the analytic value to <1e-6 relative.
    #[test]
    fn resting_point_normal_force_balances_weight() {
        let m = 2.0;
        let k_n = 5.0e5;
        let params = ContactParams::new(k_n, 50.0, 0.0, 1e-3);
        let plane = Plane::ground(0.0);
        // Place the point exactly at the analytic equilibrium penetration and
        // verify the spring force there equals the weight, at zero velocity.
        let d_eq = m * G / k_n;
        let p = Vector2::new(0.0, -d_eq); // d_eq below the surface
        let f = contact_force(p, Vector2::zeros(), &plane, &params);
        let weight = m * G;
        // At rest (ḋ = 0) the force is purely the spring: kₙ·d = m·g.
        assert!(f.x.abs() < 1e-9, "no tangential force at rest, got {}", f.x);
        let rel = (f.y - weight).abs() / weight;
        assert!(
            rel < 1e-6,
            "normal force {} vs weight {} (rel {rel:e})",
            f.y,
            weight
        );
        // And the equilibrium penetration itself: solve kₙ·d = m·g → d = m·g/kₙ.
        let d_from_force = f.y / k_n;
        assert!(
            (d_from_force - d_eq).abs() / d_eq < 1e-6,
            "penetration {d_from_force} vs analytic {d_eq}"
        );
    }

    /// A free body dropped onto the ground and integrated to rest settles at
    /// the analytic steady penetration `m·g/kₙ` (closes the loop dynamically,
    /// not just at the analytic point).
    #[test]
    fn dropped_body_settles_at_analytic_penetration() {
        let m = 1.5;
        let k_n = 2.0e4;
        let params = ContactParams::new(k_n, 80.0, 0.0, 1e-3);
        let plane = Plane::ground(0.0);
        // A single body, contact point at its CoM (so no torque).
        let mut body = Body::new(m, 0.05, Vector2::new(0.0, 0.02)); // start above ground
        let contacts = [BodyContact {
            local: Vector2::zeros(),
            plane,
            params,
        }];
        let dt = 1.0e-4;
        for _ in 0..200_000 {
            let state = ContactBodyState::from(&body);
            let w = body_contact_wrench(&state, &contacts);
            // Newton: a = (gravity + contact)/m ; symplectic Euler.
            let ay = -G + w.force.y / m;
            body.vel.y += ay * dt;
            body.pos.y += body.vel.y * dt;
        }
        assert!(body.vel.y.abs() < 1e-3, "did not settle, v {}", body.vel.y);
        let d = plane.penetration(body.pos);
        let d_eq = m * G / k_n;
        assert!(
            (d - d_eq).abs() / d_eq < 1e-3,
            "settled penetration {d} vs analytic {d_eq}"
        );
    }

    /// (2) A point clear of the plane (no penetration) produces **exactly
    /// zero** contact force — including a point moving fast toward the plane
    /// but not yet touching.
    #[test]
    fn no_penetration_gives_zero_force() {
        let params = ContactParams::new(1.0e5, 100.0, 0.7, 1e-3);
        let plane = Plane::ground(0.0);
        // Above the ground, moving down hard — still no force until it crosses.
        let above = Vector2::new(0.3, 0.5);
        let f = contact_force(above, Vector2::new(0.0, -10.0), &plane, &params);
        assert_eq!(f, Vector2::zeros(), "force above ground must be zero");
        // Exactly on the surface (depth 0) is also not penetrating.
        let on = Vector2::new(0.0, 0.0);
        let f0 = contact_force(on, Vector2::zeros(), &plane, &params);
        assert_eq!(f0, Vector2::zeros(), "force on surface must be zero");
    }

    /// The normal damper is **unilateral**: a point penetrating only shallowly
    /// but separating fast must not be *pulled* back — the force clamps to zero
    /// rather than going negative (no glue).
    #[test]
    fn separating_contact_never_pulls() {
        let params = ContactParams::new(1.0e3, 500.0, 0.0, 1e-3);
        let plane = Plane::ground(0.0);
        let p = Vector2::new(0.0, -1e-4); // barely inside
                                          // Moving up (separating) fast: cₙ·ḋ is large and negative, would make
                                          // kₙ·d + cₙ·ḋ < 0. Must clamp to zero, not pull down.
        let f = contact_force(p, Vector2::new(0.0, 5.0), &plane, &params);
        assert_eq!(
            f,
            Vector2::zeros(),
            "separating contact must not pull, {f:?}"
        );
    }

    /// (3a) A tangential push **below** the friction limit `μ·fₙ` keeps a
    /// resting block **stuck** — its tangential velocity stays near zero.
    #[test]
    fn tangential_push_below_limit_sticks() {
        let m = 1.0;
        let k_n = 1.0e5;
        let mu = 0.5;
        let params = ContactParams::new(k_n, 100.0, mu, 1e-4);
        let plane = Plane::ground(0.0);
        let mut body = Body::new(m, 0.05, Vector2::new(0.0, -m * G / k_n)); // resting
        let contacts = [BodyContact {
            local: Vector2::zeros(),
            plane,
            params,
        }];
        // Horizontal applied force well under μ·m·g (the slip threshold).
        let f_applied = 0.4 * mu * m * G;
        let dt = 1.0e-4;
        for _ in 0..50_000 {
            let state = ContactBodyState::from(&body);
            let w = body_contact_wrench(&state, &contacts);
            let ax = (f_applied + w.force.x) / m;
            let ay = (-G * m + w.force.y) / m;
            body.vel.x += ax * dt;
            body.vel.y += ay * dt;
            body.pos.x += body.vel.x * dt;
            body.pos.y += body.vel.y * dt;
        }
        // Stuck: it should barely creep (regularized friction allows a tiny
        // viscous crawl, but no runaway sliding).
        assert!(
            body.vel.x.abs() < 0.02,
            "below-limit push should stick, vx {}",
            body.vel.x
        );
        assert!(body.pos.x.abs() < 0.01, "crept too far, x {}", body.pos.x);
    }

    /// (3b) A tangential push **above** the friction limit makes the block
    /// **slide**, and once sliding the kinetic friction force equals `μ·fₙ`
    /// opposing the motion (matched analytically).
    #[test]
    fn tangential_push_above_limit_slides_with_kinetic_friction() {
        let m = 1.0;
        let k_n = 1.0e5;
        let mu = 0.5;
        let v_eps = 1e-4;
        let params = ContactParams::new(k_n, 100.0, mu, v_eps);
        let plane = Plane::ground(0.0);
        let mut body = Body::new(m, 0.05, Vector2::new(0.0, -m * G / k_n));
        let contacts = [BodyContact {
            local: Vector2::zeros(),
            plane,
            params,
        }];
        // Push at 3× the slip threshold → net along-ground accel ≈ (3−1)μg.
        let f_applied = 3.0 * mu * m * G;
        let dt = 1.0e-4;
        for _ in 0..50_000 {
            let state = ContactBodyState::from(&body);
            let w = body_contact_wrench(&state, &contacts);
            let ax = (f_applied + w.force.x) / m;
            let ay = (-G * m + w.force.y) / m;
            body.vel.x += ax * dt;
            body.vel.y += ay * dt;
            body.pos.x += body.vel.x * dt;
            body.pos.y += body.vel.y * dt;
        }
        // It is sliding fast (|v_t| ≫ v_eps), so friction has saturated.
        assert!(
            body.vel.x > 1.0,
            "above-limit push should slide, vx {}",
            body.vel.x
        );
        // Kinetic friction magnitude == μ·fₙ, opposing motion (force.x < 0).
        let state = ContactBodyState::from(&body);
        let w = body_contact_wrench(&state, &contacts);
        let f_n = k_n * plane.penetration(body.pos); // at near-zero vertical speed
        let expected_friction = mu * f_n;
        assert!(
            w.force.x < 0.0,
            "friction must oppose motion, {}",
            w.force.x
        );
        let rel = (w.force.x.abs() - expected_friction).abs() / expected_friction;
        assert!(
            rel < 1e-2,
            "kinetic friction {} vs μ·fₙ {} (rel {rel:e})",
            w.force.x.abs(),
            expected_friction
        );
    }

    /// Friction is bounded by the cone: for *any* sliding speed the tangential
    /// force magnitude never exceeds `μ·fₙ` (a property test of the `tanh`
    /// saturation).
    #[test]
    fn friction_never_exceeds_cone() {
        let params = ContactParams::new(1.0e4, 0.0, 0.6, 1e-3);
        let plane = Plane::ground(0.0);
        let p = Vector2::new(0.0, -1e-3); // fixed penetration → fixed fₙ
        let f_n = params.normal_stiffness * plane.penetration(p);
        let cone = params.friction * f_n;
        for &v in &[1e-6, 1e-3, 0.01, 0.1, 1.0, 10.0, 100.0] {
            let f = contact_force(p, Vector2::new(v, 0.0), &plane, &params);
            // Tangential component is f.x here (normal is +y).
            assert!(
                f.x.abs() <= cone + 1e-9,
                "friction {} exceeded cone {} at v={v}",
                f.x.abs(),
                cone
            );
            // And it opposes the motion.
            assert!(f.x <= 0.0, "friction must oppose +x motion at v={v}");
        }
    }

    /// (4) A body **dropped** onto the ground with damping settles and **never
    /// gains energy**: the peak total mechanical energy after release never
    /// exceeds the release energy (penalty damping only dissipates).
    #[test]
    fn dropped_body_never_gains_energy() {
        let m = 1.0;
        let k_n = 1.0e4;
        let c_n = 60.0;
        let params = ContactParams::new(k_n, c_n, 0.0, 1e-3);
        let plane = Plane::ground(0.0);
        let h0 = 0.05;
        let mut body = Body::new(m, 0.05, Vector2::new(0.0, h0));
        let contacts = [BodyContact {
            local: Vector2::zeros(),
            plane,
            params,
        }];
        // Total mechanical energy = KE + gravity PE + contact spring PE.
        // Spring PE while penetrating depth d is ½·kₙ·d² (a conservative store);
        // including it makes "never gains energy" a clean statement — the damper
        // is the only non-conservative term and it can only remove energy.
        let energy = |b: &Body| -> f64 {
            let ke = 0.5 * m * b.vel.norm_squared();
            let pe_g = m * G * b.pos.y;
            let d = plane.penetration(b.pos);
            let pe_c = if d > 0.0 { 0.5 * k_n * d * d } else { 0.0 };
            ke + pe_g + pe_c
        };
        let e0 = energy(&body);
        let dt = 1.0e-4;
        let mut max_e = e0;
        for _ in 0..300_000 {
            let state = ContactBodyState::from(&body);
            let w = body_contact_wrench(&state, &contacts);
            let ay = -G + w.force.y / m;
            body.vel.y += ay * dt;
            body.pos.y += body.vel.y * dt;
            max_e = max_e.max(energy(&body));
        }
        // Never gains energy (small symplectic/explicit tolerance).
        assert!(
            max_e <= e0 + 1e-6 * e0.abs().max(1.0),
            "energy grew: peak {max_e} vs release {e0}"
        );
        // And it has actually settled near the analytic resting penetration.
        assert!(body.vel.y.abs() < 1e-2, "did not settle, v {}", body.vel.y);
        let d_eq = m * G / k_n;
        assert!(
            (plane.penetration(body.pos) - d_eq).abs() / d_eq < 5e-2,
            "did not settle near analytic penetration"
        );
    }

    /// A contact point off the centre of mass produces a **torque** about the
    /// CoM (the wrench reduction carries the moment arm `r × f`), and a point
    /// on the CoM produces none.
    #[test]
    fn off_center_contact_produces_torque() {
        let params = ContactParams::new(1.0e4, 10.0, 0.0, 1e-3);
        let plane = Plane::ground(0.0);
        // Body CoM above ground; a contact point offset +x and pushed into the
        // ground (place the corner below the surface).
        let state = ContactBodyState {
            pos: Vector2::new(0.0, 0.01),
            angle: 0.0,
            vel: Vector2::zeros(),
            omega: 0.0,
        };
        let corner = Vector2::new(0.5, -0.02); // local point, 0.02 below ground in world
        let contacts = [BodyContact {
            local: corner,
            plane,
            params,
        }];
        let w = body_contact_wrench(&state, &contacts);
        assert!(w.force.y > 0.0, "should push up, {}", w.force.y);
        // Force is +y at arm +x → torque r × f = arm.x·f.y > 0 (CCW).
        assert!(
            w.torque > 0.0,
            "off-centre normal force should torque, {}",
            w.torque
        );
        // Same contact on the CoM → no torque.
        let centred = [BodyContact {
            local: Vector2::new(0.0, -0.02),
            plane,
            params,
        }];
        let w0 = body_contact_wrench(&state, &centred);
        assert!(
            w0.torque.abs() < 1e-9,
            "centred contact torque {}",
            w0.torque
        );
    }

    /// `From<&Body>` produces the same state as building it by hand (the bridge
    /// to the existing [`crate::Body`] is faithful).
    #[test]
    fn body_state_bridge_is_faithful() {
        let mut b = Body::new(2.0, 0.3, Vector2::new(1.0, 2.0));
        b.angle = 0.5;
        b.vel = Vector2::new(-0.1, 0.2);
        b.omega = 1.3;
        let s = ContactBodyState::from(&b);
        assert_eq!(s.pos, b.pos);
        assert_eq!(s.angle, b.angle);
        assert_eq!(s.vel, b.vel);
        assert_eq!(s.omega, b.omega);
    }

    /// A non-flat (tilted) plane still gives the right penetration sign and a
    /// force along its normal.
    #[test]
    fn tilted_plane_force_is_along_normal() {
        // Plane through origin, normal pointing up-and-right (45°), normalised
        // by `new`.
        let plane = Plane::new(Vector2::zeros(), Vector2::new(1.0, 1.0));
        let params = ContactParams::new(1.0e4, 0.0, 0.0, 1e-3);
        // A point on the inside (opposite the normal) penetrates.
        let p = Vector2::new(-0.1, -0.1);
        let d = plane.penetration(p);
        assert!(d > 0.0, "should penetrate, depth {d}");
        let f = contact_force(p, Vector2::zeros(), &plane, &params);
        // Force is purely normal (frictionless) → parallel to (1,1)/√2.
        let n = plane.normal;
        let along = f.dot(&n);
        let perp_comp = (f - n * along).norm();
        assert!(perp_comp < 1e-9, "force off-normal by {perp_comp}");
        assert!(along > 0.0, "force should push out along +n");
    }

    // ---- fail-loud configuration tests ----

    #[test]
    #[should_panic(expected = "normal_stiffness")]
    fn rejects_nonpositive_stiffness() {
        ContactParams::new(0.0, 1.0, 0.5, 1e-3);
    }

    #[test]
    #[should_panic(expected = "normal_damping")]
    fn rejects_negative_damping() {
        ContactParams::new(1.0e4, -1.0, 0.5, 1e-3);
    }

    #[test]
    #[should_panic(expected = "friction")]
    fn rejects_negative_friction() {
        ContactParams::new(1.0e4, 1.0, -0.1, 1e-3);
    }

    #[test]
    #[should_panic(expected = "friction_vel_scale")]
    fn rejects_nonpositive_vel_scale() {
        ContactParams::new(1.0e4, 1.0, 0.5, 0.0);
    }

    #[test]
    #[should_panic(expected = "normal_stiffness")]
    fn rejects_nan_stiffness() {
        ContactParams::new(f64::NAN, 1.0, 0.5, 1e-3);
    }

    #[test]
    #[should_panic(expected = "normal")]
    fn rejects_zero_plane_normal() {
        Plane::new(Vector2::zeros(), Vector2::zeros());
    }

    #[test]
    #[should_panic(expected = "normal")]
    fn rejects_nan_plane_normal() {
        Plane::new(Vector2::zeros(), Vector2::new(f64::NAN, 1.0));
    }
}
