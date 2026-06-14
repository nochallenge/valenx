//! # valenx-mbd
//!
//! A native **planar (2-D) multibody-dynamics** solver — the time-domain
//! dynamics engine valenx was missing (its assembly model is kinematics-only,
//! and dynamics otherwise route to the MuJoCo adapter). Build a mechanism from
//! **rigid bodies** (mass + rotational inertia), **force elements** (gravity,
//! linear springs, dampers) and **holonomic constraints** (revolute pins,
//! prismatic sliders), and [`System::step`] advances it through time.
//!
//! ## The method
//!
//! This is the genuine constrained-multibody formulation, the same one a
//! production code (Adams) is built on:
//!
//! - Generalised coordinates `q = (x, y, θ)` per body; a block-diagonal mass
//!   matrix `M = diag(m, m, I)`.
//! - Applied forces assemble into a generalised force `Q = (Fx, Fy, τ)`.
//! - Each holonomic constraint `C(q) = 0` contributes its Jacobian `Cq`; the
//!   **constrained equations of motion** are the saddle-point (KKT) system
//!   `[M Cqᵀ; Cq 0]·[q̈; λ] = [Q; γ]`, with `λ` the Lagrange multipliers
//!   (constraint forces) and `γ` the acceleration-level right-hand side
//!   (centripetal term + **Baumgarte stabilization** `−2α·Ċ − β²·C` to stop
//!   constraint drift).
//! - Time integration is **semi-implicit (symplectic) Euler** — velocities
//!   first, then positions — which keeps the energy of conservative systems
//!   bounded.
//!
//! ## Validated
//!
//! The tests check it against closed-form results: a physical pendulum's small
//! oscillation period `2π√(I_pivot/(m g d))`, a spring–mass natural frequency
//! `√(k/m)`, spring–mass **energy conservation**, and a double pendulum's
//! bounded energy.
//!
//! ## Honest scope
//!
//! Planar (2-D), rigid bodies only, with revolute pins, prismatic sliders,
//! springs and dampers — a real v1, **research / preliminary-design grade**. It
//! is a step toward, not an equal of, a production multibody code (Adams): no
//! 3-D spatial dynamics (quaternions / 3×3 inertia), no contact, friction,
//! flexible bodies, cylindrical/gear joints, or an implicit stiff integrator
//! yet. The prismatic slider here is the ground variant (a body on a fixed
//! world axis); a body-to-body prismatic and the rest are the documented next
//! steps.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use nalgebra::{DMatrix, DVector, Vector2};

/// Baumgarte position/velocity stabilization gains.
const BAUMGARTE_ALPHA: f64 = 20.0;
const BAUMGARTE_BETA: f64 = 20.0;

/// A planar rigid body — generalised coordinates `(x, y, θ)` and their rates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Body {
    /// Mass (kg).
    pub mass: f64,
    /// Rotational inertia about the centre of mass (kg·m²).
    pub inertia: f64,
    /// Centre-of-mass position (m).
    pub pos: Vector2<f64>,
    /// Orientation (rad).
    pub angle: f64,
    /// Centre-of-mass velocity (m/s).
    pub vel: Vector2<f64>,
    /// Angular velocity (rad/s).
    pub omega: f64,
}

impl Body {
    /// A body at rest at `pos`, angle 0.
    pub fn new(mass: f64, inertia: f64, pos: Vector2<f64>) -> Self {
        Self {
            mass,
            inertia,
            pos,
            angle: 0.0,
            vel: Vector2::zeros(),
            omega: 0.0,
        }
    }
}

/// Where a force/constraint attaches: a body-fixed local point, or a fixed
/// point in the world.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Anchor {
    /// A point fixed in body `index`'s local frame.
    Body {
        /// Body index.
        index: usize,
        /// Attachment point in the body's local frame (m).
        local: Vector2<f64>,
    },
    /// A point fixed in the world (m).
    World(Vector2<f64>),
}

/// A force element applied to the system.
#[derive(Debug, Clone, Copy)]
pub enum Force {
    /// Uniform gravitational acceleration applied to every body's centre of
    /// mass.
    Gravity(Vector2<f64>),
    /// A linear spring + damper between two anchors: axial force
    /// `k·(len − rest) + c·(rate of length change)`.
    Spring {
        /// First attachment.
        a: Anchor,
        /// Second attachment.
        b: Anchor,
        /// Stiffness `k` (N/m).
        stiffness: f64,
        /// Natural (unstretched) length (m).
        rest_length: f64,
        /// Damping `c` (N·s/m).
        damping: f64,
    },
}

/// A holonomic constraint.
#[derive(Debug, Clone, Copy)]
pub enum Constraint {
    /// Revolute pin coupling two anchors so they coincide (2 scalar
    /// constraints). A `World`–`Body` pair pins a body to the ground; a
    /// `Body`–`Body` pair is a revolute joint.
    Pin(Anchor, Anchor),
    /// Prismatic (slider) joint pinning a body to a fixed **world axis**: it may
    /// only translate along that axis, with its perpendicular offset and its
    /// orientation both locked (2 scalar constraints). This is the ground
    /// variant — the classic block-on-a-ramp slider.
    Prismatic {
        /// Index of the sliding body.
        body: usize,
        /// A point the slide axis passes through, in world coordinates (m).
        axis_point: Vector2<f64>,
        /// Slide-axis direction in the world; need not be unit (normalised
        /// internally). Must be non-zero.
        axis: Vector2<f64>,
        /// The body orientation the joint holds fixed (rad) — a slider does not
        /// rotate, so this is normally the body's initial angle.
        ref_angle: f64,
    },
}

/// A multibody system: bodies, force elements, constraints and the clock.
#[derive(Debug, Clone, Default)]
pub struct System {
    /// The rigid bodies.
    pub bodies: Vec<Body>,
    /// The force elements.
    pub forces: Vec<Force>,
    /// The constraints.
    pub constraints: Vec<Constraint>,
    /// Simulation time (s).
    pub time: f64,
}

fn rot(angle: f64, v: Vector2<f64>) -> Vector2<f64> {
    let (s, c) = angle.sin_cos();
    Vector2::new(c * v.x - s * v.y, s * v.x + c * v.y)
}

fn perp(v: Vector2<f64>) -> Vector2<f64> {
    Vector2::new(-v.y, v.x)
}

fn cross2(a: Vector2<f64>, b: Vector2<f64>) -> f64 {
    a.x * b.y - a.y * b.x
}

impl System {
    /// World position of an anchor, and (for a body anchor) the rotated local
    /// arm `rp = R(θ)·local`.
    fn anchor_world(&self, a: Anchor) -> (Vector2<f64>, Vector2<f64>, Option<usize>) {
        match a {
            Anchor::World(p) => (p, Vector2::zeros(), None),
            Anchor::Body { index, local } => {
                let b = &self.bodies[index];
                let rp = rot(b.angle, local);
                (b.pos + rp, rp, Some(index))
            }
        }
    }

    /// Velocity of an anchor point.
    fn anchor_velocity(&self, a: Anchor) -> Vector2<f64> {
        match a {
            Anchor::World(_) => Vector2::zeros(),
            Anchor::Body { index, local } => {
                let b = &self.bodies[index];
                let rp = rot(b.angle, local);
                b.vel + b.omega * perp(rp)
            }
        }
    }

    /// Assemble the generalised applied force `Q` (length `3n`).
    fn applied_forces(&self) -> DVector<f64> {
        let n = self.bodies.len();
        let mut q = DVector::zeros(3 * n);
        for f in &self.forces {
            match *f {
                Force::Gravity(g) => {
                    for (i, b) in self.bodies.iter().enumerate() {
                        q[3 * i] += b.mass * g.x;
                        q[3 * i + 1] += b.mass * g.y;
                    }
                }
                Force::Spring {
                    a,
                    b,
                    stiffness,
                    rest_length,
                    damping,
                } => {
                    let (pa, rpa, ia) = self.anchor_world(a);
                    let (pb, rpb, ib) = self.anchor_world(b);
                    let d = pb - pa;
                    let len = d.norm();
                    if len < 1e-12 {
                        continue;
                    }
                    let dir = d / len;
                    let rate = (self.anchor_velocity(b) - self.anchor_velocity(a)).dot(&dir);
                    let mag = stiffness * (len - rest_length) + damping * rate;
                    let force_on_a = dir * mag; // pulls a toward b when stretched
                    if let Some(i) = ia {
                        q[3 * i] += force_on_a.x;
                        q[3 * i + 1] += force_on_a.y;
                        q[3 * i + 2] += cross2(rpa, force_on_a);
                    }
                    if let Some(i) = ib {
                        q[3 * i] -= force_on_a.x;
                        q[3 * i + 1] -= force_on_a.y;
                        q[3 * i + 2] -= cross2(rpb, force_on_a);
                    }
                }
            }
        }
        q
    }

    /// Number of scalar constraint rows (each pin or prismatic joint is 2).
    fn constraint_rows(&self) -> usize {
        self.constraints
            .iter()
            .map(|c| match c {
                Constraint::Pin(..) => 2,
                Constraint::Prismatic { .. } => 2,
            })
            .sum()
    }

    /// Assemble the constraint Jacobian `Cq`, value `C` and RHS `γ`.
    fn assemble_constraints(&self) -> (DMatrix<f64>, DVector<f64>) {
        let n = self.bodies.len();
        let m = self.constraint_rows();
        let mut cq = DMatrix::zeros(m, 3 * n);
        let mut gamma = DVector::zeros(m);
        let mut row = 0;
        for con in &self.constraints {
            match *con {
                Constraint::Pin(a, b) => {
                    let (pa, rpa, ia) = self.anchor_world(a);
                    let (pb, rpb, ib) = self.anchor_world(b);
                    let c = pa - pb; // = 0 when satisfied
                    let cdot = self.anchor_velocity(a) - self.anchor_velocity(b);
                    // Jacobian: ∂(pa−pb)/∂q.
                    if let Some(i) = ia {
                        cq[(row, 3 * i)] += 1.0;
                        cq[(row, 3 * i + 2)] += -rpa.y;
                        cq[(row + 1, 3 * i + 1)] += 1.0;
                        cq[(row + 1, 3 * i + 2)] += rpa.x;
                    }
                    if let Some(i) = ib {
                        cq[(row, 3 * i)] += -1.0;
                        cq[(row, 3 * i + 2)] += rpb.y;
                        cq[(row + 1, 3 * i + 1)] += -1.0;
                        cq[(row + 1, 3 * i + 2)] += -rpb.x;
                    }
                    // γ = (ω_a²·rp_a − ω_b²·rp_b) − 2α·Ċ − β²·C (centripetal + Baumgarte).
                    let mut centripetal = Vector2::zeros();
                    if let Anchor::Body { index, .. } = a {
                        centripetal += self.bodies[index].omega.powi(2) * rpa;
                    }
                    if let Anchor::Body { index, .. } = b {
                        centripetal -= self.bodies[index].omega.powi(2) * rpb;
                    }
                    let g = centripetal - 2.0 * BAUMGARTE_ALPHA * cdot - BAUMGARTE_BETA.powi(2) * c;
                    gamma[row] = g.x;
                    gamma[row + 1] = g.y;
                    row += 2;
                }
                Constraint::Prismatic {
                    body,
                    axis_point,
                    axis,
                    ref_angle,
                } => {
                    let bd = &self.bodies[body];
                    // Unit slide direction `u` and its world-fixed normal `nrm`.
                    let u = axis.normalize();
                    let nrm = perp(u);
                    // Row 0: no perpendicular drift — C = nrm·(pos − axis_point).
                    // Both nrm and the constraint are linear in q, so the
                    // acceleration-level RHS is pure Baumgarte (no centripetal
                    // term): γ = −2α·Ċ − β²·C, with Ċ = nrm·vel.
                    let c_perp = nrm.dot(&(bd.pos - axis_point));
                    let cdot_perp = nrm.dot(&bd.vel);
                    cq[(row, 3 * body)] += nrm.x;
                    cq[(row, 3 * body + 1)] += nrm.y;
                    gamma[row] =
                        -2.0 * BAUMGARTE_ALPHA * cdot_perp - BAUMGARTE_BETA.powi(2) * c_perp;
                    // Row 1: no rotation — C = θ − ref_angle, Ċ = ω.
                    let c_ang = bd.angle - ref_angle;
                    cq[(row + 1, 3 * body + 2)] += 1.0;
                    gamma[row + 1] =
                        -2.0 * BAUMGARTE_ALPHA * bd.omega - BAUMGARTE_BETA.powi(2) * c_ang;
                    row += 2;
                }
            }
        }
        (cq, gamma)
    }

    /// Solve for the generalised accelerations (length `3n`).
    fn accelerations(&self) -> DVector<f64> {
        let n = self.bodies.len();
        let q = self.applied_forces();
        let m = self.constraint_rows();
        if m == 0 {
            // Unconstrained: M is diagonal, so divide directly.
            let mut a = DVector::zeros(3 * n);
            for (i, b) in self.bodies.iter().enumerate() {
                a[3 * i] = q[3 * i] / b.mass;
                a[3 * i + 1] = q[3 * i + 1] / b.mass;
                a[3 * i + 2] = q[3 * i + 2] / b.inertia;
            }
            return a;
        }
        let (cq, gamma) = self.assemble_constraints();
        let size = 3 * n + m;
        let mut kkt = DMatrix::zeros(size, size);
        for (i, b) in self.bodies.iter().enumerate() {
            kkt[(3 * i, 3 * i)] = b.mass;
            kkt[(3 * i + 1, 3 * i + 1)] = b.mass;
            kkt[(3 * i + 2, 3 * i + 2)] = b.inertia;
        }
        for r in 0..m {
            for col in 0..3 * n {
                let v = cq[(r, col)];
                if v != 0.0 {
                    kkt[(3 * n + r, col)] = v; // Cq
                    kkt[(col, 3 * n + r)] = v; // Cqᵀ
                }
            }
        }
        let mut rhs = DVector::zeros(size);
        rhs.rows_mut(0, 3 * n).copy_from(&q);
        rhs.rows_mut(3 * n, m).copy_from(&gamma);
        let sol = kkt.lu().solve(&rhs).unwrap_or_else(|| DVector::zeros(size));
        sol.rows(0, 3 * n).into_owned()
    }

    /// Advance the system by one time step `dt` (semi-implicit Euler).
    pub fn step(&mut self, dt: f64) {
        let a = self.accelerations();
        for (i, b) in self.bodies.iter_mut().enumerate() {
            b.vel.x += a[3 * i] * dt;
            b.vel.y += a[3 * i + 1] * dt;
            b.omega += a[3 * i + 2] * dt;
            b.pos.x += b.vel.x * dt;
            b.pos.y += b.vel.y * dt;
            b.angle += b.omega * dt;
        }
        self.time += dt;
    }

    /// Total mechanical energy (kinetic + gravitational + spring potential), J.
    pub fn energy(&self) -> f64 {
        let mut e = 0.0;
        for b in &self.bodies {
            e += 0.5 * b.mass * b.vel.norm_squared() + 0.5 * b.inertia * b.omega * b.omega;
        }
        for f in &self.forces {
            match *f {
                Force::Gravity(g) => {
                    for b in &self.bodies {
                        e += -b.mass * g.dot(&b.pos);
                    }
                }
                Force::Spring {
                    a,
                    b,
                    stiffness,
                    rest_length,
                    ..
                } => {
                    let len = (self.anchor_world(a).0 - self.anchor_world(b).0).norm();
                    e += 0.5 * stiffness * (len - rest_length).powi(2);
                }
            }
        }
        e
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Run `steps` of size `dt`, returning the angle (or x) sampler results.
    fn period_from_zero_crossings(samples: &[(f64, f64)]) -> Option<f64> {
        let mut crossings = Vec::new();
        for w in samples.windows(2) {
            let (t0, y0) = w[0];
            let (t1, y1) = w[1];
            if y0 == 0.0 || (y0 < 0.0) != (y1 < 0.0) {
                // Linear-interpolate the crossing time.
                let frac = if (y1 - y0).abs() > 1e-15 {
                    -y0 / (y1 - y0)
                } else {
                    0.0
                };
                crossings.push(t0 + frac * (t1 - t0));
            }
        }
        // One full period = every other zero crossing.
        (crossings.len() >= 3).then(|| crossings[2] - crossings[0])
    }

    #[test]
    fn physical_pendulum_matches_analytic_period() {
        let (m, i_cg, d, theta0) = (1.0, 0.10, 1.0, 0.05);
        let g = 9.81;
        // Pin a body (CG at distance d below the pivot) to the world origin.
        let local_pin = Vector2::new(0.0, d); // pivot is `d` above the CG
        let pos = -rot(theta0, local_pin); // keep the pin on the world origin
        let mut body = Body::new(m, i_cg, pos);
        body.angle = theta0;
        let mut sys = System {
            bodies: vec![body],
            forces: vec![Force::Gravity(Vector2::new(0.0, -g))],
            constraints: vec![Constraint::Pin(
                Anchor::Body {
                    index: 0,
                    local: local_pin,
                },
                Anchor::World(Vector2::zeros()),
            )],
            time: 0.0,
        };
        let dt = 5.0e-4;
        let mut samples = Vec::new();
        for _ in 0..12_000 {
            samples.push((sys.time, sys.bodies[0].angle));
            sys.step(dt);
        }
        let measured = period_from_zero_crossings(&samples).expect("oscillation");
        let i_pivot = i_cg + m * d * d;
        let analytic = 2.0 * PI * (i_pivot / (m * g * d)).sqrt();
        assert!(
            (measured - analytic).abs() / analytic < 0.02,
            "period {measured} vs analytic {analytic}"
        );
    }

    #[test]
    fn spring_mass_matches_natural_frequency() {
        let (m, k, rest) = (1.0, 10.0, 1.0);
        // Body displaced 0.1 m along x from the rest position; spring to ground.
        let mut sys = System {
            bodies: vec![Body::new(m, 0.1, Vector2::new(rest + 0.1, 0.0))],
            forces: vec![Force::Spring {
                a: Anchor::World(Vector2::zeros()),
                b: Anchor::Body {
                    index: 0,
                    local: Vector2::zeros(),
                },
                stiffness: k,
                rest_length: rest,
                damping: 0.0,
            }],
            constraints: vec![],
            time: 0.0,
        };
        let dt = 2.0e-4;
        let mut samples = Vec::new();
        for _ in 0..30_000 {
            samples.push((sys.time, sys.bodies[0].pos.x - rest));
            sys.step(dt);
        }
        let measured = period_from_zero_crossings(&samples).expect("oscillation");
        let analytic = 2.0 * PI * (m / k).sqrt();
        assert!(
            (measured - analytic).abs() / analytic < 0.02,
            "period {measured} vs analytic {analytic}"
        );
    }

    #[test]
    fn undamped_spring_mass_conserves_energy() {
        let mut sys = System {
            bodies: vec![Body::new(1.0, 0.1, Vector2::new(1.2, 0.0))],
            forces: vec![Force::Spring {
                a: Anchor::World(Vector2::zeros()),
                b: Anchor::Body {
                    index: 0,
                    local: Vector2::zeros(),
                },
                stiffness: 12.0,
                rest_length: 1.0,
                damping: 0.0,
            }],
            constraints: vec![],
            time: 0.0,
        };
        let e0 = sys.energy();
        let mut max_dev = 0.0_f64;
        for _ in 0..40_000 {
            sys.step(1.0e-4);
            max_dev = max_dev.max((sys.energy() - e0).abs() / e0);
        }
        assert!(max_dev < 0.01, "symplectic energy drift {max_dev}");
    }

    #[test]
    fn double_pendulum_energy_stays_bounded() {
        // Two hanging links (body0 pinned to ground at the top, body1 pinned to
        // the bottom of body0), released from a small consistent tilt at rest —
        // a genuine body-to-body revolute chain. A horizontal release would be
        // a violent max-energy start that an explicit integrator can't hold.
        let pivot = Vector2::new(0.0, 0.0);
        let (t0, t1) = (0.12_f64, 0.12_f64); // small tilts from straight down
        let cg0 = pivot + rot(t0, Vector2::new(0.0, -0.5));
        let end0 = pivot + rot(t0, Vector2::new(0.0, -1.0));
        let cg1 = end0 + rot(t1, Vector2::new(0.0, -0.5));
        let mut b0 = Body::new(1.0, 0.05, cg0);
        b0.angle = t0;
        let mut b1 = Body::new(1.0, 0.05, cg1);
        b1.angle = t1;
        let mut sys = System {
            bodies: vec![b0, b1],
            forces: vec![Force::Gravity(Vector2::new(0.0, -9.81))],
            constraints: vec![
                Constraint::Pin(
                    Anchor::Body {
                        index: 0,
                        local: Vector2::new(0.0, 0.5),
                    },
                    Anchor::World(pivot),
                ),
                Constraint::Pin(
                    Anchor::Body {
                        index: 0,
                        local: Vector2::new(0.0, -0.5),
                    },
                    Anchor::Body {
                        index: 1,
                        local: Vector2::new(0.0, 0.5),
                    },
                ),
            ],
            time: 0.0,
        };
        let e0 = sys.energy();
        let mut max_dev = 0.0_f64;
        for _ in 0..20_000 {
            sys.step(1.0e-4);
            max_dev = max_dev.max((sys.energy() - e0).abs() / e0.abs());
        }
        // Chaotic but Hamiltonian: energy must stay bounded (Baumgarte adds a
        // little non-conservation, hence a looser tolerance than the
        // unconstrained symplectic case).
        assert!(max_dev < 0.05, "double-pendulum energy drift {max_dev}");
    }

    #[test]
    fn prismatic_incline_accelerates_at_g_sin_theta() {
        // A frictionless block on a ramp inclined θ above horizontal must
        // accelerate down-slope at exactly g·sinθ — the textbook result.
        let g = 9.81;
        let theta = 30.0_f64.to_radians();
        let u = Vector2::new(theta.cos(), theta.sin()); // up-slope direction
        let mut sys = System {
            bodies: vec![Body::new(2.0, 0.1, Vector2::zeros())],
            forces: vec![Force::Gravity(Vector2::new(0.0, -g))],
            constraints: vec![Constraint::Prismatic {
                body: 0,
                axis_point: Vector2::zeros(),
                axis: u,
                ref_angle: 0.0,
            }],
            time: 0.0,
        };
        let dt = 1.0e-4;
        let n = 5_000;
        for _ in 0..n {
            sys.step(dt);
        }
        let b = &sys.bodies[0];
        let t = n as f64 * dt;
        let nrm = perp(u);
        // Gravity pulls down-slope (−u), so the along-axis velocity is negative.
        let v_along = b.vel.dot(&u);
        let v_perp = b.vel.dot(&nrm);
        let expected = -g * theta.sin() * t;
        assert!(
            (v_along - expected).abs() < 0.01 * expected.abs(),
            "along-axis speed {v_along} vs g·sinθ·t {expected}"
        );
        assert!(v_perp.abs() < 1e-3, "perpendicular drift speed {v_perp}");
        assert!(
            b.omega.abs() < 1e-6,
            "slider should not rotate, ω {}",
            b.omega
        );
    }

    #[test]
    fn prismatic_horizontal_slider_stays_put_under_gravity() {
        // A horizontal slider: gravity is fully carried by the joint normal
        // (g·sin0 = 0 along the axis), so the body must not move.
        let g = 9.81;
        let mut sys = System {
            bodies: vec![Body::new(1.0, 0.1, Vector2::new(0.5, 0.0))],
            forces: vec![Force::Gravity(Vector2::new(0.0, -g))],
            constraints: vec![Constraint::Prismatic {
                body: 0,
                axis_point: Vector2::zeros(),
                axis: Vector2::new(1.0, 0.0),
                ref_angle: 0.0,
            }],
            time: 0.0,
        };
        let start = sys.bodies[0].pos;
        for _ in 0..5_000 {
            sys.step(1.0e-4);
        }
        let b = &sys.bodies[0];
        assert!(
            (b.pos - start).norm() < 1e-3,
            "slider drifted {:?}",
            b.pos - start
        );
        assert!(b.vel.norm() < 1e-3, "slider gained speed {:?}", b.vel);
    }

    #[test]
    fn prismatic_slider_coasts_at_constant_velocity() {
        // No force along the axis → an initial glide is preserved (momentum).
        let mut body = Body::new(1.0, 0.1, Vector2::zeros());
        body.vel = Vector2::new(0.3, 0.0);
        let mut sys = System {
            bodies: vec![body],
            forces: vec![],
            constraints: vec![Constraint::Prismatic {
                body: 0,
                axis_point: Vector2::zeros(),
                axis: Vector2::new(1.0, 0.0),
                ref_angle: 0.0,
            }],
            time: 0.0,
        };
        for _ in 0..1_000 {
            sys.step(1.0e-3); // 1.0 s total
        }
        let b = &sys.bodies[0];
        assert!(
            (b.vel.x - 0.3).abs() < 1e-6,
            "axial speed changed: {}",
            b.vel.x
        );
        assert!(
            (b.pos.x - 0.3).abs() < 1e-3,
            "x after 1 s: {} (want 0.3)",
            b.pos.x
        );
        assert!(
            b.pos.y.abs() < 1e-6 && b.vel.y.abs() < 1e-6,
            "left the axis"
        );
    }

    #[test]
    fn is_deterministic() {
        let make = || System {
            bodies: vec![Body::new(1.0, 0.1, Vector2::new(1.1, 0.0))],
            forces: vec![Force::Spring {
                a: Anchor::World(Vector2::zeros()),
                b: Anchor::Body {
                    index: 0,
                    local: Vector2::zeros(),
                },
                stiffness: 10.0,
                rest_length: 1.0,
                damping: 0.1,
            }],
            constraints: vec![],
            time: 0.0,
        };
        let (mut a, mut b) = (make(), make());
        for _ in 0..5_000 {
            a.step(1.0e-4);
            b.step(1.0e-4);
        }
        assert_eq!(a.bodies[0].pos, b.bodies[0].pos);
    }
}
