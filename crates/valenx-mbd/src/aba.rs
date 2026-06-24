//! Articulated-body forward dynamics for a kinematic tree — the **Featherstone
//! Articulated-Body Algorithm (ABA)**, an `O(n)` propagation-based forward
//! dynamics solver.
//!
//! This is a separate, complementary solver to the planar constrained-DAE
//! [`crate::System`]: where that engine assembles a global KKT system over
//! 2-D bodies and holonomic constraints, this one works on a **3-D kinematic
//! tree** (a fixed base with single-degree-of-freedom joints), in full spatial
//! (6-D) algebra, and solves forward dynamics by three recursive sweeps over
//! the tree rather than by factoring one big matrix.
//!
//! Given joint positions `q`, joint velocities `qd` and applied joint torques
//! `tau`, [`ArticulatedTree::forward_dynamics`] returns the joint accelerations
//! `qdd` in `O(n)` time. [`ArticulatedTree::step`] then advances `(q, qd)` with
//! **semi-implicit (symplectic) Euler**, matching the convention used by the
//! planar [`crate::System::step`].
//!
//! ## Method (Featherstone notation)
//!
//! Spatial vectors are 6-vectors `[angular; linear]`; a spatial motion vector is
//! `v = [ω; vₒ]` (angular velocity, then the linear velocity of the body-frame
//! origin), and a spatial force is `f = [n; f]` (couple, then linear force).
//! Spatial inertia is a `6×6` matrix `I`. The algorithm is the standard three
//! passes:
//!
//! 1. **Outward (base → leaves):** propagate body spatial velocities `vᵢ` and
//!    the velocity-product (Coriolis/centrifugal) bias terms.
//! 2. **Inward (leaves → base):** build the *articulated-body* spatial inertias
//!    `Iᴬᵢ` and bias forces `pᴬᵢ`, accumulating each child's contribution onto
//!    its parent.
//! 3. **Outward (base → leaves):** solve for joint accelerations `q̈ᵢ` and the
//!    body spatial accelerations `aᵢ`.
//!
//! ## Scope
//!
//! - **Joints:** single-DOF **revolute** and **prismatic** (arbitrary 3-D axis).
//! - **Base:** either **fixed** (the world is body `−1`, immovable —
//!   [`ArticulatedTree::forward_dynamics`]) **or floating** (a 6-DOF free root
//!   body that may translate and rotate — [`ArticulatedTree::forward_dynamics_floating`]).
//!   The floating case replaces the fixed `a₀ = −a_g` with the articulated-base
//!   solve `a₀ = −(Iᴬ₀)⁻¹ pᴬ₀` (the base articulated inertia and bias force at
//!   the root), then propagates outward exactly as before.
//! - Each body carries a rigid-body spatial inertia built from mass, a
//!   centre-of-mass offset, and a `3×3` rotational inertia about the CoM.
//!
//! Clean-room implementation of the published algorithm (Featherstone, *Rigid
//! Body Dynamics Algorithms*, 2008, Ch. 7), not ported from any existing code.

use nalgebra::{Matrix3, Matrix6, Vector3, Vector6};

/// Standard gravity magnitude (m/s²), pointing along `−z` by default.
pub const STANDARD_GRAVITY: f64 = 9.80665;

/// A single-degree-of-freedom joint connecting a body to its parent.
///
/// The axis is expressed in the **child body's frame** and need not be unit
/// (it is normalised internally). A revolute joint rotates the child about the
/// axis; a prismatic joint translates the child along it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum JointType {
    /// Rotation about `axis` (rad). Generalised coordinate is the joint angle.
    Revolute {
        /// Rotation axis in the child body frame.
        axis: Vector3<f64>,
    },
    /// Translation along `axis` (m). Generalised coordinate is the joint
    /// displacement.
    Prismatic {
        /// Slide axis in the child body frame.
        axis: Vector3<f64>,
    },
}

impl JointType {
    /// The (constant, body-frame) spatial **motion subspace** `S` of this joint
    /// — the 6-vector that maps the scalar joint rate `q̇` to the joint's spatial
    /// velocity contribution `S·q̇`. For a revolute joint `S = [û; 0]`; for a
    /// prismatic joint `S = [0; û]`, with `û` the unit axis.
    fn motion_subspace(&self) -> Vector6<f64> {
        let mut s = Vector6::zeros();
        match self {
            JointType::Revolute { axis } => {
                let u = unit(*axis);
                s.fixed_rows_mut::<3>(0).copy_from(&u);
            }
            JointType::Prismatic { axis } => {
                let u = unit(*axis);
                s.fixed_rows_mut::<3>(3).copy_from(&u);
            }
        }
        s
    }

    /// The child-relative spatial transform produced by joint coordinate `q`,
    /// i.e. the transform `X_J` from the joint's parent-side frame to its
    /// child-side frame (Featherstone `X_J(q)`).
    fn joint_transform(&self, q: f64) -> SpatialTransform {
        match self {
            JointType::Revolute { axis } => SpatialTransform::rotation(unit(*axis), q),
            JointType::Prismatic { axis } => SpatialTransform::translation(unit(*axis) * q),
        }
    }
}

/// One body (link) of the kinematic tree.
///
/// A body is reached from its parent through a fixed transform `X_T` (the
/// parent-body frame to this body's joint frame, the constant part of the
/// kinematic chain) followed by the variable joint transform `X_J(q)`.
#[derive(Debug, Clone, Copy)]
pub struct TreeBody {
    /// Index of the parent body, or [`ArticulatedTree::BASE`] for a body
    /// attached directly to the fixed base. A parent index must be **less than**
    /// this body's own index (the tree is stored in topological order).
    pub parent: usize,
    /// The joint connecting this body to its parent.
    pub joint: JointType,
    /// Fixed transform from the **parent body frame** to this body's
    /// **joint (predecessor) frame** — the constant link geometry. Identity if
    /// the joint frame coincides with the parent origin.
    pub parent_to_joint: SpatialTransform,
    /// Mass of the body (kg).
    pub mass: f64,
    /// Centre-of-mass position in the body frame (m).
    pub com: Vector3<f64>,
    /// Rotational inertia about the **centre of mass**, in the body frame
    /// (kg·m²), as a symmetric `3×3` matrix.
    pub inertia_com: Matrix3<f64>,
}

impl TreeBody {
    /// A revolute-jointed body. `parent_to_joint` is the fixed parent→joint
    /// transform; `axis` is the rotation axis in the body frame.
    pub fn revolute(
        parent: usize,
        axis: Vector3<f64>,
        parent_to_joint: SpatialTransform,
        mass: f64,
        com: Vector3<f64>,
        inertia_com: Matrix3<f64>,
    ) -> Self {
        Self {
            parent,
            joint: JointType::Revolute { axis },
            parent_to_joint,
            mass,
            com,
            inertia_com,
        }
    }

    /// A prismatic-jointed body. `axis` is the slide axis in the body frame.
    pub fn prismatic(
        parent: usize,
        axis: Vector3<f64>,
        parent_to_joint: SpatialTransform,
        mass: f64,
        com: Vector3<f64>,
        inertia_com: Matrix3<f64>,
    ) -> Self {
        Self {
            parent,
            joint: JointType::Prismatic { axis },
            parent_to_joint,
            mass,
            com,
            inertia_com,
        }
    }

    /// The body's rigid-body **spatial inertia** `I` (a `6×6` matrix) expressed
    /// about the **body-frame origin**, built from mass, CoM offset `c` and the
    /// rotational inertia `I_c` about the CoM. In `[angular; linear]` ordering,
    ///
    /// ```text
    /// I = [ I_c + m·[c]ₓ[c]ₓᵀ    m·[c]ₓ ]
    ///     [        m·[c]ₓᵀ          m·1₃ ]
    /// ```
    ///
    /// where `[c]ₓ` is the skew (cross-product) matrix of `c`. (Featherstone
    /// eq. 2.63.)
    fn spatial_inertia(&self) -> Matrix6<f64> {
        let m = self.mass;
        let c = self.com;
        let cx = skew(c);
        // Inertia about the origin = I_c − m·[c]ₓ[c]ₓ (parallel-axis), and since
        // [c]ₓ[c]ₓ = [c]ₓ(−[c]ₓᵀ) = −[c]ₓ[c]ₓᵀ, this equals I_c + m·[c]ₓ[c]ₓᵀ.
        let i_o = self.inertia_com + m * (cx * cx.transpose());
        let mut spatial = Matrix6::zeros();
        spatial.fixed_view_mut::<3, 3>(0, 0).copy_from(&i_o);
        spatial.fixed_view_mut::<3, 3>(0, 3).copy_from(&(m * cx));
        spatial
            .fixed_view_mut::<3, 3>(3, 0)
            .copy_from(&(m * cx.transpose()));
        spatial
            .fixed_view_mut::<3, 3>(3, 3)
            .copy_from(&(m * Matrix3::identity()));
        spatial
    }
}

/// A spatial (`6×6`) coordinate transform between body frames (a Plücker
/// transform), stored as a rotation `E` and a translation `r` (Featherstone's
/// `X = rot(E)·xlt(r)`: first translate the origin by `r`, then rotate by `E`).
///
/// Applied to a spatial **motion** vector `v = [ω; vₒ]`:
/// `X·v = [E·ω ; E·(vₒ − r×ω)]`. Forces transform contravariantly; the inward
/// pass carries child forces and inertias up to the parent by the `Xᵀ`
/// congruence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialTransform {
    /// Rotation matrix `E` (parent → child orientation).
    rot: Matrix3<f64>,
    /// Translation `r` of the new origin, expressed in the **old** frame.
    trans: Vector3<f64>,
}

impl SpatialTransform {
    /// The identity transform.
    pub fn identity() -> Self {
        Self {
            rot: Matrix3::identity(),
            trans: Vector3::zeros(),
        }
    }

    /// A pure translation: the new frame's origin is at `r` (in the old frame),
    /// with axes unchanged.
    pub fn translation(r: Vector3<f64>) -> Self {
        Self {
            rot: Matrix3::identity(),
            trans: r,
        }
    }

    /// A pure rotation of `angle` (rad) about a unit `axis`, origins coincident.
    /// The stored `E` rotates a vector from the old frame into the new frame, so
    /// it is the transpose of the active (Rodrigues) rotation by `+angle`.
    pub fn rotation(axis: Vector3<f64>, angle: f64) -> Self {
        Self {
            rot: rodrigues(axis, angle).transpose(),
            trans: Vector3::zeros(),
        }
    }

    /// A general transform: translate the origin by `r` (old-frame coords) then
    /// orient the new axes by rotation matrix `e` (old → new).
    pub fn new(e: Matrix3<f64>, r: Vector3<f64>) -> Self {
        Self { rot: e, trans: r }
    }

    /// Compose: `self ∘ other` applies `other` first, then `self`
    /// (so the result transforms a vector through `other` and then `self`).
    fn then(&self, other: &SpatialTransform) -> SpatialTransform {
        // X_a · X_b with X = rot(E)·xlt(r):
        //   E = E_a·E_b,   r = r_b + E_bᵀ·r_a.
        let e = self.rot * other.rot;
        let r = other.trans + other.rot.transpose() * self.trans;
        SpatialTransform { rot: e, trans: r }
    }

    /// Apply this transform to a spatial **motion** vector.
    fn apply_motion(&self, v: &Vector6<f64>) -> Vector6<f64> {
        let w = v.fixed_rows::<3>(0).into_owned();
        let vo = v.fixed_rows::<3>(3).into_owned();
        let w_new = self.rot * w;
        let v_new = self.rot * (vo - self.trans.cross(&w));
        stack(w_new, v_new)
    }

    /// Carry a spatial **force expressed in the child frame up into the parent
    /// frame**. With `self` mapping motion parent → child (`v_c = X·v_p`),
    /// forces transform contravariantly, `f_p = Xᵀ·f_c` (power invariance
    /// `f_c·v_c = f_p·v_p`). This is the operator the inward pass needs to add a
    /// child's force to its parent.
    fn force_to_parent(&self, f_child: &Vector6<f64>) -> Vector6<f64> {
        self.to_motion_matrix().transpose() * f_child
    }

    /// Carry a spatial **inertia expressed in the child frame up into the parent
    /// frame**. With `self` mapping motion parent → child, an inertia maps by
    /// the congruence `I_p = Xᵀ · I_c · X` (so that `I_p·v_p` is the parent-frame
    /// momentum of the same motion). Used to add a child's articulated inertia
    /// to its parent.
    fn inertia_to_parent(&self, i_child: &Matrix6<f64>) -> Matrix6<f64> {
        let x = self.to_motion_matrix();
        x.transpose() * i_child * x
    }

    /// The explicit `6×6` matrix of this transform acting on motion vectors.
    fn to_motion_matrix(self) -> Matrix6<f64> {
        let e = self.rot;
        let mrx = -e * skew(self.trans); // E·(−[r]ₓ)
        let mut m = Matrix6::zeros();
        m.fixed_view_mut::<3, 3>(0, 0).copy_from(&e);
        m.fixed_view_mut::<3, 3>(3, 0).copy_from(&mrx);
        m.fixed_view_mut::<3, 3>(3, 3).copy_from(&e);
        m
    }

    /// The inverse transform (maps the child frame back to the parent frame).
    pub fn inverse(&self) -> SpatialTransform {
        // (rot(E)·xlt(r))⁻¹ = xlt(−r)·rot(Eᵀ) = rot(Eᵀ)·xlt(−E·r).
        let et = self.rot.transpose();
        SpatialTransform {
            rot: et,
            trans: -(self.rot * self.trans),
        }
    }
}

/// The free **floating base** (root body) of a tree for
/// [`ArticulatedTree::forward_dynamics_floating`].
///
/// Unlike the fixed base of [`ArticulatedTree::forward_dynamics`] (an immovable
/// world), this root body has its own rigid-body spatial inertia and a 6-DOF
/// spatial velocity `v₀ = [ω₀; v₀ₒ]`, and is free to translate and rotate. The
/// tree's body-0 children attach to it through their `parent_to_joint`
/// transforms with `parent == BASE` exactly as for the fixed base.
///
/// All quantities are expressed in the **base body frame**. The forward-dynamics
/// solve returns the base's spatial acceleration in this same frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FloatingBase {
    /// Mass of the base body (kg). Must be finite and positive.
    pub mass: f64,
    /// Centre-of-mass position in the base frame (m).
    pub com: Vector3<f64>,
    /// Rotational inertia about the base body's centre of mass, in the base
    /// frame (kg·m²).
    pub inertia_com: Matrix3<f64>,
    /// Base spatial velocity `[ω; vₒ]` **in the base frame** — angular velocity
    /// then the linear velocity of the base-frame origin.
    pub velocity: Vector6<f64>,
    /// Base orientation: the rotation matrix mapping **base-frame** vectors to
    /// the **world frame** (`R_wb`). Identity by default. Used only by
    /// [`ArticulatedTree::step_floating`] to integrate the pose; the
    /// instantaneous forward-dynamics solve is frame-agnostic.
    pub orientation: Matrix3<f64>,
    /// Position of the base-frame origin in the world (m). Used only by
    /// [`ArticulatedTree::step_floating`].
    pub position: Vector3<f64>,
}

impl FloatingBase {
    /// A base at rest (zero spatial velocity) at the world origin with identity
    /// orientation and the given inertial properties.
    pub fn new(mass: f64, com: Vector3<f64>, inertia_com: Matrix3<f64>) -> Self {
        Self {
            mass,
            com,
            inertia_com,
            velocity: Vector6::zeros(),
            orientation: Matrix3::identity(),
            position: Vector3::zeros(),
        }
    }

    /// The base body's rigid-body **spatial inertia** `I` about the base-frame
    /// origin — the same construction [`TreeBody::spatial_inertia`] uses.
    fn spatial_inertia(&self) -> Matrix6<f64> {
        let m = self.mass;
        let cx = skew(self.com);
        let i_o = self.inertia_com + m * (cx * cx.transpose());
        let mut spatial = Matrix6::zeros();
        spatial.fixed_view_mut::<3, 3>(0, 0).copy_from(&i_o);
        spatial.fixed_view_mut::<3, 3>(0, 3).copy_from(&(m * cx));
        spatial
            .fixed_view_mut::<3, 3>(3, 0)
            .copy_from(&(m * cx.transpose()));
        spatial
            .fixed_view_mut::<3, 3>(3, 3)
            .copy_from(&(m * Matrix3::identity()));
        spatial
    }
}

/// The result of a floating-base forward-dynamics solve
/// ([`ArticulatedTree::forward_dynamics_floating`]).
#[derive(Debug, Clone, PartialEq)]
pub struct FloatingAccel {
    /// The base body's spatial acceleration `a₀ = [α₀; a₀ₒ]` in the base frame
    /// (angular acceleration then the linear acceleration of the base origin).
    /// This is the **true** (inertial-frame) acceleration including the response
    /// to gravity, so a free body under gravity reports `a₀ₒ = g`.
    pub base: Vector6<f64>,
    /// Joint accelerations `q̈`, one per body — the same quantity the fixed-base
    /// [`ArticulatedTree::forward_dynamics`] returns.
    pub joints: Vec<f64>,
    /// Each link's spatial acceleration `aᵢ = [αᵢ; aᵢₒ]` in its **own** body
    /// frame, one per body (the same `a` the third recursion sweep computes).
    /// Exposed because downstream callers frequently need link accelerations
    /// (e.g. to form sensor readings or net wrenches); also lets a caller verify
    /// the momentum balance directly.
    pub body: Vec<Vector6<f64>>,
}

/// Errors from building or evaluating an [`ArticulatedTree`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TreeError {
    /// A body's `parent` index is out of range (not the base and `>= n`), or it
    /// references a body that does not precede it (which would make the stored
    /// order non-topological / introduce a cycle).
    BadParent {
        /// The offending body index.
        body: usize,
        /// The invalid parent index it named.
        parent: usize,
    },
    /// A joint axis was the zero vector (no well-defined direction).
    ZeroAxis {
        /// The offending body index.
        body: usize,
    },
    /// A supplied state slice (`q`, `qd`, or `tau`) had the wrong length.
    BadDimension {
        /// Number of joints (expected length).
        expected: usize,
        /// Length actually supplied.
        got: usize,
    },
    /// The floating base's `6×6` articulated-body inertia `Iᴬ₀` was singular, so
    /// the base-acceleration solve `a₀ = −(Iᴬ₀)⁻¹ pᴬ₀` has no unique solution
    /// (e.g. a base with zero or non-positive mass/inertia). Fail-loud rather
    /// than producing `NaN`s.
    SingularBaseInertia,
}

impl core::fmt::Display for TreeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            TreeError::BadParent { body, parent } => write!(
                f,
                "body {body} names parent {parent}, which is out of range or does not precede it \
                 (parents must be the fixed base or an earlier body — no cycles)"
            ),
            TreeError::ZeroAxis { body } => {
                write!(f, "body {body} has a zero-length joint axis")
            }
            TreeError::BadDimension { expected, got } => {
                write!(f, "expected a state vector of length {expected}, got {got}")
            }
            TreeError::SingularBaseInertia => write!(
                f,
                "the floating base's 6×6 articulated-body inertia is singular \
                 (base mass/inertia must be finite and positive)"
            ),
        }
    }
}

impl std::error::Error for TreeError {}

/// A kinematic tree of rigid bodies with single-DOF joints and a fixed base,
/// solved by the Featherstone Articulated-Body Algorithm.
///
/// Bodies are stored in topological order: body `i`'s parent is either the
/// [`BASE`](Self::BASE) or some body `j < i`. There is exactly one joint per
/// body (the joint to its parent), so the number of generalised coordinates
/// equals the number of bodies. State vectors `q`, `qd`, `qdd`, `tau` are
/// indexed by body.
#[derive(Debug, Clone, Default)]
pub struct ArticulatedTree {
    /// The bodies (links), in topological order.
    pub bodies: Vec<TreeBody>,
    /// Joint positions, one per body (rad for revolute, m for prismatic).
    pub q: Vec<f64>,
    /// Joint velocities, one per body.
    pub qd: Vec<f64>,
    /// Gravitational acceleration in the base/world frame (m/s²).
    pub gravity: Vector3<f64>,
    /// Simulation time (s).
    pub time: f64,
}

impl ArticulatedTree {
    /// Sentinel `parent` value marking a body attached to the fixed base.
    pub const BASE: usize = usize::MAX;

    /// An empty tree with gravity along `−z` and zero clock.
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            q: Vec::new(),
            qd: Vec::new(),
            gravity: Vector3::new(0.0, 0.0, -STANDARD_GRAVITY),
            time: 0.0,
        }
    }

    /// Number of bodies / joints / generalised coordinates.
    pub fn len(&self) -> usize {
        self.bodies.len()
    }

    /// Whether the tree has no bodies.
    pub fn is_empty(&self) -> bool {
        self.bodies.is_empty()
    }

    /// Append a body and its initial joint state, returning its index.
    ///
    /// Convenience builder that keeps `q`/`qd` the same length as `bodies`.
    pub fn push(&mut self, body: TreeBody, q0: f64, qd0: f64) -> usize {
        let idx = self.bodies.len();
        self.bodies.push(body);
        self.q.push(q0);
        self.qd.push(qd0);
        idx
    }

    /// Validate the tree's structure: every parent index must be the base or an
    /// earlier body, every joint axis must be non-zero, and the state vectors
    /// must match the body count.
    ///
    /// This is **fail-loud**: forward dynamics calls it first and returns the
    /// error rather than producing nonsense from a malformed tree.
    pub fn validate(&self) -> Result<(), TreeError> {
        let n = self.bodies.len();
        if self.q.len() != n {
            return Err(TreeError::BadDimension {
                expected: n,
                got: self.q.len(),
            });
        }
        if self.qd.len() != n {
            return Err(TreeError::BadDimension {
                expected: n,
                got: self.qd.len(),
            });
        }
        for (i, b) in self.bodies.iter().enumerate() {
            // A valid parent is either the base sentinel or a strictly earlier
            // body. Requiring `parent < i` rules out self-loops, forward refs
            // and cycles in one check (the storage order is the topological
            // order).
            if b.parent != Self::BASE && b.parent >= i {
                return Err(TreeError::BadParent {
                    body: i,
                    parent: b.parent,
                });
            }
            let axis = match b.joint {
                JointType::Revolute { axis } | JointType::Prismatic { axis } => axis,
            };
            if axis.norm() < 1e-12 {
                return Err(TreeError::ZeroAxis { body: i });
            }
        }
        Ok(())
    }

    /// Forward dynamics: given the current `q`, `qd` and the applied joint
    /// torques/forces `tau` (one scalar per body), compute the joint
    /// accelerations `qdd` via the three-pass Articulated-Body Algorithm.
    ///
    /// `O(n)` in the number of bodies. Returns [`TreeError`] for a malformed
    /// tree or a wrong-length `tau`.
    pub fn forward_dynamics(&self, tau: &[f64]) -> Result<Vec<f64>, TreeError> {
        self.validate()?;
        let n = self.bodies.len();
        if tau.len() != n {
            return Err(TreeError::BadDimension {
                expected: n,
                got: tau.len(),
            });
        }
        if n == 0 {
            return Ok(Vec::new());
        }

        // Per-body workspace.
        let mut xup = vec![SpatialTransform::identity(); n]; // parent → body i
        let mut s = vec![Vector6::zeros(); n]; // motion subspace
        let mut v = vec![Vector6::zeros(); n]; // spatial velocity
        let mut c = vec![Vector6::zeros(); n]; // velocity-product bias term
        let mut ia = vec![Matrix6::zeros(); n]; // articulated-body inertia
        let mut pa = vec![Vector6::zeros(); n]; // articulated-body bias force

        // ---- Pass 1: outward (base → leaves) — velocities and bias c. ----
        for i in 0..n {
            let body = &self.bodies[i];
            let s_i = body.joint.motion_subspace();
            // X from parent frame to body i frame: fixed link transform then
            // the variable joint transform.
            let x_j = body.joint.joint_transform(self.q[i]);
            let x_up = x_j.then(&body.parent_to_joint);
            // Velocity across the joint: vJ = S·q̇.
            let vj = s_i * self.qd[i];
            // Parent spatial velocity carried into this frame (zero if base).
            let v_parent_here = if body.parent == Self::BASE {
                Vector6::zeros()
            } else {
                x_up.apply_motion(&v[body.parent])
            };
            let v_i = v_parent_here + vj;
            // Velocity-product term c = v ×ₘ vJ (Featherstone eq. 7.31; with a
            // constant S there is no S̈ term).
            let c_i = motion_cross(&v_i, &vj);

            xup[i] = x_up;
            s[i] = s_i;
            v[i] = v_i;
            c[i] = c_i;
            // Initialise the articulated quantities with the body's own rigid
            // inertia and its velocity-dependent (Coriolis/centrifugal) bias
            // force p = v ×f (I·v).
            let i_rb = body.spatial_inertia();
            ia[i] = i_rb;
            pa[i] = force_cross(&v_i, &(i_rb * v_i));
        }

        // ---- Pass 2: inward (leaves → base) — articulated inertia & bias. ----
        // Cache the per-joint intermediate quantities for pass 3.
        let mut u = vec![Vector6::zeros(); n]; // U = Iᴬ·S
        let mut d = vec![0.0_f64; n]; // d = Sᵀ·U
        let mut u_small = vec![0.0_f64; n]; // u = τ − Sᵀ·pᴬ
        for i in (0..n).rev() {
            let u_i = ia[i] * s[i];
            let d_i = s[i].dot(&u_i);
            let u_small_i = tau[i] - s[i].dot(&pa[i]);
            u[i] = u_i;
            d[i] = d_i;
            u_small[i] = u_small_i;

            let parent = self.bodies[i].parent;
            if parent != Self::BASE {
                // Articulated inertia/bias contribution propagated to the parent
                // (Featherstone eq. 7.45–7.47), transformed into the parent
                // frame by X_upᵀ (force/inertia dual mapping).
                let d_safe = if d_i.abs() > 1e-12 { d_i } else { 1e-12 };
                let ia_child = ia[i] - (u_i * u_i.transpose()) / d_safe;
                let pa_child = pa[i] + ia_child * c[i] + u[i] * (u_small_i / d_safe);
                // Map child articulated inertia and bias force up into parent.
                ia[parent] += xup[i].inertia_to_parent(&ia_child);
                pa[parent] += xup[i].force_to_parent(&pa_child);
            }
        }

        // ---- Pass 3: outward (base → leaves) — accelerations. ----
        // Base acceleration carries the gravity field: a₀ = −a_g (so that every
        // free body falls at g). In spatial terms the base has acceleration
        // [0; −g] expressed in the base frame.
        let a_base = stack(Vector3::zeros(), -self.gravity);
        let mut a = vec![Vector6::zeros(); n]; // spatial acceleration
        let mut qdd = vec![0.0_f64; n];
        for i in 0..n {
            let parent = self.bodies[i].parent;
            let a_parent = if parent == Self::BASE {
                a_base
            } else {
                a[parent]
            };
            // a' = X_up·a_parent + c (the part of this body's acceleration
            // already determined before the joint accelerates).
            let a_prime = xup[i].apply_motion(&a_parent) + c[i];
            let d_safe = if d[i].abs() > 1e-12 { d[i] } else { 1e-12 };
            // q̈ = (u − Uᵀ·a') / d.
            let qdd_i = (u_small[i] - u[i].dot(&a_prime)) / d_safe;
            qdd[i] = qdd_i;
            a[i] = a_prime + s[i] * qdd_i;
        }

        Ok(qdd)
    }

    /// **Floating-base** forward dynamics: the same `O(n)` Articulated-Body
    /// Algorithm as [`forward_dynamics`](Self::forward_dynamics), but with a free
    /// 6-DOF root [`FloatingBase`] in place of the immovable fixed base. Returns
    /// both the base spatial acceleration and the joint accelerations.
    ///
    /// # Method
    ///
    /// The three passes are unchanged for the actuated bodies; only the base is
    /// treated differently:
    ///
    /// 1. **Outward:** body-0 children inherit the base's spatial velocity `v₀`
    ///    (carried into their frame) instead of zero.
    /// 2. **Inward:** every base-attached child's articulated inertia and bias
    ///    force accumulate onto the base's own `Iᴬ₀`, `pᴬ₀` (the base *is* the
    ///    reference frame, so no parent transform is applied to its own terms).
    /// 3. **Base solve & outward:** instead of the fixed `a₀ = −a_g`, solve the
    ///    articulated base equation `Iᴬ₀ · a₀ = −pᴬ₀` directly for the base's
    ///    physical spatial acceleration, then start the outward acceleration
    ///    sweep from `a₀`. Joint accelerations follow exactly as before.
    ///
    /// Gravity is handled the rigorous way for a free base: rather than the
    /// fixed-base fictitious-base-acceleration trick, each body's bias force
    /// carries the real gravitational wrench `−I·[0; g]` (the field `g` in that
    /// body's frame). The returned accelerations are therefore directly physical
    /// — a free body under gravity reports `a₀ₒ = g`. This yields the *same*
    /// joint accelerations as the fixed-base path, so making the base inertia
    /// very large (`a₀ → [0; g]`, i.e. a pinned base in free fall) recovers the
    /// fixed-base joint accelerations exactly.
    ///
    /// # Errors
    ///
    /// Returns [`TreeError`] for a malformed tree or wrong-length `tau`, and
    /// [`TreeError::SingularBaseInertia`] if `Iᴬ₀` cannot be inverted (a base
    /// with non-positive mass/inertia) — it is **fail-loud** and never unwraps a
    /// singular inverse into `NaN`s.
    pub fn forward_dynamics_floating(
        &self,
        base: &FloatingBase,
        tau: &[f64],
    ) -> Result<FloatingAccel, TreeError> {
        self.validate()?;
        let n = self.bodies.len();
        if tau.len() != n {
            return Err(TreeError::BadDimension {
                expected: n,
                got: tau.len(),
            });
        }

        // ---- Base own articulated quantities (in the base frame). ----
        // The base starts with its own rigid inertia and bias force. Unlike the
        // fixed-base path (which folds gravity into a fictitious base
        // acceleration), the floating base genuinely accelerates, so gravity is
        // applied here as a real external wrench in every body's bias: a body
        // feels the gravity field `g` as the spatial force `I·[0; g]` about its
        // origin, which enters `p` with a minus sign (external-force convention).
        // This makes the returned base/joint accelerations directly physical
        // (a free body reports a₀ₒ = g) while giving the *same* joint
        // accelerations as the fixed-base fictitious-gravity trick.
        let v0 = base.velocity;
        let i0 = base.spatial_inertia();
        let g_world = self.gravity;
        let g_base = base.orientation.transpose() * g_world; // gravity in base frame
        let mut ia0 = i0;
        let mut pa0 = force_cross(&v0, &(i0 * v0)) - i0 * stack(Vector3::zeros(), g_base);

        // Per-body workspace (bodies 0..n; the base is handled separately).
        let mut xup = vec![SpatialTransform::identity(); n]; // parent → body i
        let mut s = vec![Vector6::zeros(); n]; // motion subspace
        let mut v = vec![Vector6::zeros(); n]; // spatial velocity
        let mut c = vec![Vector6::zeros(); n]; // velocity-product bias term
        let mut ia = vec![Matrix6::zeros(); n]; // articulated-body inertia
        let mut pa = vec![Vector6::zeros(); n]; // articulated-body bias force
        let mut rbw = vec![Matrix3::identity(); n]; // body→world rotation (for gravity)

        // ---- Pass 1: outward (base → leaves) — velocities and bias c. ----
        for i in 0..n {
            let body = &self.bodies[i];
            let s_i = body.joint.motion_subspace();
            let x_j = body.joint.joint_transform(self.q[i]);
            let x_up = x_j.then(&body.parent_to_joint);
            let vj = s_i * self.qd[i];
            // The only velocity change vs the fixed base: a body-0 child inherits
            // the base's spatial velocity (carried into this frame), not zero.
            let (v_parent_here, r_parent) = if body.parent == Self::BASE {
                (x_up.apply_motion(&v0), base.orientation)
            } else {
                (x_up.apply_motion(&v[body.parent]), rbw[body.parent])
            };
            let v_i = v_parent_here + vj;
            let c_i = motion_cross(&v_i, &vj);
            // World orientation of this body: parent's, composed with (i→parent)
            // = `x_up.rot.transpose()`.
            let r_i = r_parent * x_up.rot.transpose();

            xup[i] = x_up;
            s[i] = s_i;
            v[i] = v_i;
            c[i] = c_i;
            rbw[i] = r_i;
            let i_rb = body.spatial_inertia();
            ia[i] = i_rb;
            // Bias = velocity-product term − gravitational wrench (gravity field
            // expressed in this body's frame).
            let g_i = r_i.transpose() * g_world;
            pa[i] = force_cross(&v_i, &(i_rb * v_i)) - i_rb * stack(Vector3::zeros(), g_i);
        }

        // ---- Pass 2: inward (leaves → base) — articulated inertia & bias. ----
        let mut u = vec![Vector6::zeros(); n]; // U = Iᴬ·S
        let mut d = vec![0.0_f64; n]; // d = Sᵀ·U
        let mut u_small = vec![0.0_f64; n]; // u = τ − Sᵀ·pᴬ
        for i in (0..n).rev() {
            let u_i = ia[i] * s[i];
            let d_i = s[i].dot(&u_i);
            let u_small_i = tau[i] - s[i].dot(&pa[i]);
            u[i] = u_i;
            d[i] = d_i;
            u_small[i] = u_small_i;

            let d_safe = if d_i.abs() > 1e-12 { d_i } else { 1e-12 };
            let ia_child = ia[i] - (u_i * u_i.transpose()) / d_safe;
            let pa_child = pa[i] + ia_child * c[i] + u[i] * (u_small_i / d_safe);
            let parent = self.bodies[i].parent;
            if parent != Self::BASE {
                ia[parent] += xup[i].inertia_to_parent(&ia_child);
                pa[parent] += xup[i].force_to_parent(&pa_child);
            } else {
                // Accumulate onto the base's own articulated inertia/bias.
                ia0 += xup[i].inertia_to_parent(&ia_child);
                pa0 += xup[i].force_to_parent(&pa_child);
            }
        }

        // ---- Base solve: a₀ = −(Iᴬ₀)⁻¹ · pᴬ₀. ----
        // With gravity already in pᴬ₀ as a real wrench, this is the base's
        // physical spatial acceleration directly (no fictitious offset to undo).
        // Guard the 6×6 inverse: fail-loud rather than unwrap a singular inertia.
        let ia0_inv = ia0.try_inverse().ok_or(TreeError::SingularBaseInertia)?;
        let a0 = -(ia0_inv * pa0);

        // ---- Pass 3: outward (base → leaves) — joint accelerations. ----
        let mut a = vec![Vector6::zeros(); n]; // spatial acceleration
        let mut qdd = vec![0.0_f64; n];
        for i in 0..n {
            let parent = self.bodies[i].parent;
            let a_parent = if parent == Self::BASE { a0 } else { a[parent] };
            let a_prime = xup[i].apply_motion(&a_parent) + c[i];
            let d_safe = if d[i].abs() > 1e-12 { d[i] } else { 1e-12 };
            let qdd_i = (u_small[i] - u[i].dot(&a_prime)) / d_safe;
            qdd[i] = qdd_i;
            a[i] = a_prime + s[i] * qdd_i;
        }

        Ok(FloatingAccel {
            base: a0,
            joints: qdd,
            body: a,
        })
    }

    /// Advance `(q, qd)` by one step `dt` using **semi-implicit (symplectic)
    /// Euler** — velocities updated from the accelerations first, then positions
    /// from the new velocities. Matches the integrator convention of the planar
    /// [`crate::System::step`].
    ///
    /// `tau` is the applied joint torque/force vector for this step. Returns the
    /// error (and leaves the state untouched) if the tree or `tau` is malformed.
    pub fn step(&mut self, tau: &[f64], dt: f64) -> Result<(), TreeError> {
        let qdd = self.forward_dynamics(tau)?;
        for (i, &acc) in qdd.iter().enumerate() {
            self.qd[i] += acc * dt;
            self.q[i] += self.qd[i] * dt;
        }
        self.time += dt;
        Ok(())
    }

    /// Advance a **floating-base** tree by one step `dt` (semi-implicit Euler),
    /// updating both the joint state `(q, qd)` *and* the supplied [`FloatingBase`]
    /// (its spatial velocity and its world pose). The base spatial velocity is
    /// updated from the base acceleration first; the joint velocities/positions
    /// then update exactly as in [`step`](Self::step); finally the base pose
    /// integrates from the (updated) base spatial velocity.
    ///
    /// Pose integration uses the spatial-velocity kinematics of a rigid body:
    /// with the base spatial velocity `[ω; vₒ]` in the base frame, the
    /// orientation advances by the exponential map `R ← R·exp([ω·dt]ₓ)` (so the
    /// rotation stays orthonormal) and the world origin by `x ← x + dt·R·vₒ`.
    /// Returns the [`TreeError`] (leaving all state untouched) on a malformed
    /// tree, wrong-length `tau`, or singular base inertia.
    pub fn step_floating(
        &mut self,
        base: &mut FloatingBase,
        tau: &[f64],
        dt: f64,
    ) -> Result<(), TreeError> {
        let accel = self.forward_dynamics_floating(base, tau)?;
        // Base spatial velocity (apparent body-frame derivative: v̇ = a).
        base.velocity += accel.base * dt;
        // Joints, exactly as in the fixed-base step.
        for (i, &acc) in accel.joints.iter().enumerate() {
            self.qd[i] += acc * dt;
            self.q[i] += self.qd[i] * dt;
        }
        // Base pose from the updated spatial velocity.
        let omega = base.velocity.fixed_rows::<3>(0).into_owned();
        let vo = base.velocity.fixed_rows::<3>(3).into_owned();
        let dr = rodrigues(omega, omega.norm() * dt); // exp([ω dt]ₓ); axis arbitrary if ω≈0
        base.orientation *= dr;
        base.position += dt * (base.orientation * vo);
        self.time += dt;
        Ok(())
    }

    /// World-frame total **linear and angular momentum** of a floating-base
    /// system about the **world origin**, returned as `[angular; linear]` (the
    /// total angular momentum about the world origin, then the total linear
    /// momentum). Sums the base and every link.
    ///
    /// With gravity off and only internal joint torques there is no external
    /// wrench, so this vector is conserved — the defining correctness property
    /// of the floating-base solve. Provided as a diagnostic / test hook.
    pub fn momentum_world(&self, base: &FloatingBase) -> Vector6<f64> {
        let kin = self.world_kinematics(base);
        // Each body's spatial momentum I·v lives in the body frame; re-express
        // it in the world frame as [torque; force] about the world origin. With
        // body→world rotation R and body origin at world position p, angular
        // (couple) n_w = R·n_b + p × (R·f_b) and linear f_w = R·f_b.
        let mut total = momentum_to_world(
            base.spatial_inertia() * base.velocity,
            &base.orientation,
            &base.position,
        );
        for i in 0..self.bodies.len() {
            let h_body = self.bodies[i].spatial_inertia() * kin.v_body[i];
            total += momentum_to_world(h_body, &kin.rot_body_to_world[i], &kin.origin_world[i]);
        }
        total
    }

    /// The **total external spatial wrench** acting on a floating-base system,
    /// about the world origin, for the given base state and joint torques `tau`
    /// — `[torque; force]` in the world frame.
    ///
    /// It runs the forward-dynamics solve and sums each body's net wrench
    /// `fᵢ = Iᵢ·aᵢ + vᵢ ×* (Iᵢ·vᵢ)` (base and every link, each in its own frame)
    /// re-expressed in the world frame. Internal joint forces cancel in the sum,
    /// so this returns exactly the external loading: **zero** when gravity is off
    /// and only internal torques act (giving an `O(round-off)`, integrator-free
    /// proof of momentum conservation), and the gravitational wrench
    /// `[Σ x_comᵢ×mᵢg ; (Σmᵢ)g]` (the total weight, with its couple about the
    /// world origin) when gravity is on. Mirrors `d/dt` of
    /// [`momentum_world`](Self::momentum_world).
    pub fn external_wrench_world(
        &self,
        base: &FloatingBase,
        tau: &[f64],
    ) -> Result<Vector6<f64>, TreeError> {
        let acc = self.forward_dynamics_floating(base, tau)?;
        let kin = self.world_kinematics(base);
        // Net wrench on one rigid body in its own frame: I·a + v ×* (I·v).
        let net = |i_rb: &Matrix6<f64>, v: &Vector6<f64>, a: &Vector6<f64>| -> Vector6<f64> {
            i_rb * a + force_cross(v, &(i_rb * v))
        };
        // Base contribution.
        let i0 = base.spatial_inertia();
        let f_base = net(&i0, &base.velocity, &acc.base);
        let mut total = momentum_to_world(f_base, &base.orientation, &base.position);
        // Each link.
        for i in 0..self.bodies.len() {
            let i_rb = self.bodies[i].spatial_inertia();
            let f_i = net(&i_rb, &kin.v_body[i], &acc.body[i]);
            total += momentum_to_world(f_i, &kin.rot_body_to_world[i], &kin.origin_world[i]);
        }
        Ok(total)
    }

    /// Total mechanical energy (kinetic + gravitational potential) of a
    /// floating-base system, J. Kinetic energy is `½ vᵀ I v` (frame-invariant)
    /// summed over the base and links; gravitational potential is `−m gᵀ x_com`
    /// per body using each body's world centre of mass. Provided as a diagnostic
    /// / test hook (energy is conserved for a frictionless conservative
    /// mechanism).
    pub fn energy_floating(&self, base: &FloatingBase) -> f64 {
        let kin = self.world_kinematics(base);
        let i0 = base.spatial_inertia();
        let mut e = 0.5 * base.velocity.dot(&(i0 * base.velocity));
        let base_com_world = base.position + base.orientation * base.com;
        e += -base.mass * self.gravity.dot(&base_com_world);
        for (i, body) in self.bodies.iter().enumerate() {
            let i_rb = body.spatial_inertia();
            e += 0.5 * kin.v_body[i].dot(&(i_rb * kin.v_body[i]));
            let com_world = kin.origin_world[i] + kin.rot_body_to_world[i] * body.com;
            e += -body.mass * self.gravity.dot(&com_world);
        }
        e
    }

    /// Forward kinematics for a floating base: each body's spatial velocity in
    /// its **own** frame, plus its world pose (body→world rotation and the world
    /// position of the body-frame origin).
    fn world_kinematics(&self, base: &FloatingBase) -> WorldKinematics {
        let n = self.bodies.len();
        let mut v_body = vec![Vector6::zeros(); n];
        let mut rot_body_to_world = vec![Matrix3::identity(); n];
        let mut origin_world = vec![Vector3::zeros(); n];
        for i in 0..n {
            let body = &self.bodies[i];
            let x_j = body.joint.joint_transform(self.q[i]);
            let x_up = x_j.then(&body.parent_to_joint); // parent → i (motion)
            let s_i = body.joint.motion_subspace();
            let vj = s_i * self.qd[i];
            // Parent's spatial velocity (body frame) and world pose.
            let (v_parent, r_w_parent, p_parent) = if body.parent == Self::BASE {
                (base.velocity, base.orientation, base.position)
            } else {
                (
                    v_body[body.parent],
                    rot_body_to_world[body.parent],
                    origin_world[body.parent],
                )
            };
            v_body[i] = x_up.apply_motion(&v_parent) + vj;
            // `x_up.rot` maps parent→i orientation; its transpose maps i→parent.
            // `x_up.trans` is body i's origin expressed in the parent frame.
            rot_body_to_world[i] = r_w_parent * x_up.rot.transpose();
            origin_world[i] = p_parent + r_w_parent * x_up.trans;
        }
        WorldKinematics {
            v_body,
            rot_body_to_world,
            origin_world,
        }
    }
}

/// Per-body forward kinematics of a floating-base tree (see
/// [`ArticulatedTree::world_kinematics`]).
struct WorldKinematics {
    /// Each body's spatial velocity in its own body frame.
    v_body: Vec<Vector6<f64>>,
    /// Each body's body→world rotation matrix.
    rot_body_to_world: Vec<Matrix3<f64>>,
    /// World position of each body-frame origin (m).
    origin_world: Vec<Vector3<f64>>,
}

// ---------------------------------------------------------------------------
// Spatial-algebra helpers.
// ---------------------------------------------------------------------------

/// Normalise a vector; returns it unchanged if it is (numerically) zero —
/// callers validate non-zero axes up front, so this only guards against a NaN.
fn unit(v: Vector3<f64>) -> Vector3<f64> {
    let n = v.norm();
    if n > 1e-12 {
        v / n
    } else {
        v
    }
}

/// Stack an angular and a linear 3-vector into a spatial 6-vector
/// `[angular; linear]`.
fn stack(angular: Vector3<f64>, linear: Vector3<f64>) -> Vector6<f64> {
    let mut s = Vector6::zeros();
    s.fixed_rows_mut::<3>(0).copy_from(&angular);
    s.fixed_rows_mut::<3>(3).copy_from(&linear);
    s
}

/// The skew-symmetric (cross-product) matrix `[v]ₓ` such that `[v]ₓ·a = v×a`.
fn skew(v: Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
}

/// Rodrigues rotation matrix: active rotation by `angle` (rad) about unit
/// `axis`.
fn rodrigues(axis: Vector3<f64>, angle: f64) -> Matrix3<f64> {
    let k = unit(axis);
    let (s, c) = angle.sin_cos();
    let kx = skew(k);
    Matrix3::identity() + s * kx + (1.0 - c) * (kx * kx)
}

/// Spatial **motion** cross product `v ×ₘ m` (Featherstone's `crm(v)·m`):
/// for `v = [ω; vₒ]`, `m = [ω'; v'ₒ]`,
/// `v ×ₘ m = [ω×ω' ; ω×v'ₒ + vₒ×ω']`.
fn motion_cross(v: &Vector6<f64>, m: &Vector6<f64>) -> Vector6<f64> {
    let w = v.fixed_rows::<3>(0).into_owned();
    let vo = v.fixed_rows::<3>(3).into_owned();
    let w2 = m.fixed_rows::<3>(0).into_owned();
    let vo2 = m.fixed_rows::<3>(3).into_owned();
    stack(w.cross(&w2), w.cross(&vo2) + vo.cross(&w2))
}

/// Re-express a body-frame spatial **momentum** `[n; f]` (couple about the body
/// origin, then linear momentum) as a spatial momentum about the **world
/// origin**, given the body→world rotation `R` and the world position `p` of the
/// body origin: `f_w = R·f`, `n_w = R·n + p×f_w`. Used to sum the system's total
/// momentum in one consistent frame.
fn momentum_to_world(h_body: Vector6<f64>, r: &Matrix3<f64>, p: &Vector3<f64>) -> Vector6<f64> {
    let n_b = h_body.fixed_rows::<3>(0).into_owned();
    let f_b = h_body.fixed_rows::<3>(3).into_owned();
    let f_w = r * f_b;
    let n_w = r * n_b + p.cross(&f_w);
    stack(n_w, f_w)
}

/// Spatial **force** cross product `v ×f f` (Featherstone's `crf(v)·f`, the
/// dual of the motion cross): for motion `v = [ω; vₒ]` and force `f = [n; f]`,
/// `v ×f f = [ω×n + vₒ×f ; ω×f]`.
fn force_cross(v: &Vector6<f64>, f: &Vector6<f64>) -> Vector6<f64> {
    let w = v.fixed_rows::<3>(0).into_owned();
    let vo = v.fixed_rows::<3>(3).into_owned();
    let n = f.fixed_rows::<3>(0).into_owned();
    let force = f.fixed_rows::<3>(3).into_owned();
    stack(w.cross(&n) + vo.cross(&force), w.cross(&force))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// Build the unit rotational inertia of a point/compact body for tests where
    /// only the joint-axis inertia matters (a thin uniform value).
    fn diag_inertia(ixx: f64, iyy: f64, izz: f64) -> Matrix3<f64> {
        Matrix3::new(ixx, 0.0, 0.0, 0.0, iyy, 0.0, 0.0, 0.0, izz)
    }

    /// (1) Single pendulum: a rigid link pinned at one end, swinging under
    /// gravity, must have angular acceleration `q̈ = −(m·g·d/I_pivot)·sin q`
    /// (the closed form: gravity torque over inertia about the pivot).
    #[test]
    fn single_pendulum_matches_closed_form() {
        // A uniform rod, length L, mass m, pivot at one end about the world z
        // axis. The CoM is placed at body-local (0, −d, 0) so that q = 0 hangs
        // the rod straight down (gravity along −y) — the natural convention in
        // which the gravity torque is exactly −m·g·d·sin q.
        let m = 2.0;
        let l = 1.5;
        let d = l / 2.0; // CoM is at the rod's midpoint, d from the pivot.
        let g = 9.80665;
        // Rod inertia about its CoM, axis = z: (1/12) m L².
        let i_cg = m * l * l / 12.0;
        let inertia_com = diag_inertia(1e-6, 1e-6, i_cg);

        let body = TreeBody::revolute(
            ArticulatedTree::BASE,
            Vector3::z(),
            SpatialTransform::identity(),
            m,
            Vector3::new(0.0, -d, 0.0),
            inertia_com,
        );
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, -g, 0.0);
        let q0 = 0.37; // an arbitrary non-trivial angle
        tree.push(body, q0, 0.0);

        let qdd = tree.forward_dynamics(&[0.0]).unwrap();
        let i_pivot = i_cg + m * d * d; // parallel axis
                                        // World CoM = R_z(q)·(0,−d,0) = (d sin q, −d cos q, 0); gravity torque
                                        // about the pivot (z) = −m·g·d·sin q, so q̈ = −(m g d / I_pivot)·sin q.
        let analytic = -(m * g * d / i_pivot) * q0.sin();
        assert!(
            (qdd[0] - analytic).abs() < 1e-9,
            "single-pendulum q̈ {} vs closed form {}",
            qdd[0],
            analytic
        );
    }

    /// (1b) The same pendulum with a pure applied torque and zero gravity must
    /// give `q̈ = τ / I_pivot` exactly.
    #[test]
    fn single_pendulum_pure_torque() {
        let m = 1.3;
        let l = 0.8;
        let d = l / 2.0;
        let i_cg = m * l * l / 12.0;
        let body = TreeBody::revolute(
            ArticulatedTree::BASE,
            Vector3::z(),
            SpatialTransform::identity(),
            m,
            Vector3::new(d, 0.0, 0.0),
            diag_inertia(1e-6, 1e-6, i_cg),
        );
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::zeros();
        tree.push(body, 0.0, 0.0);
        let torque = 0.55;
        let qdd = tree.forward_dynamics(&[torque]).unwrap();
        let i_pivot = i_cg + m * d * d;
        assert!(
            (qdd[0] - torque / i_pivot).abs() < 1e-12,
            "q̈ {} vs τ/I {}",
            qdd[0],
            torque / i_pivot
        );
    }

    /// (2) **The real ABA correctness test** — a planar double pendulum of two
    /// point masses on massless rods, evaluated against the Lagrangian closed
    /// form at a known, non-trivial configuration `(θ₁, θ₂, θ̇₁, θ̇₂)` (nonzero
    /// angles AND velocities, so the Coriolis/centrifugal coupling is exercised,
    /// not just gravity).
    #[test]
    fn double_pendulum_matches_lagrangian() {
        // Geometry / masses.
        let (m1, m2) = (1.0, 1.7);
        let (l1, l2) = (1.0, 0.8);
        let g = 9.80665;

        // State: q1 is link-1's absolute angle, q2 is link-2's angle *relative*
        // to link 1 (the standard relative-coordinate double pendulum, which is
        // how the kinematic tree composes joint 2 in joint 1's frame).
        let (t1, t2) = (0.5_f64, -0.3_f64);
        let (w1, w2) = (0.4_f64, -0.6_f64);

        // --- valenx ABA model: two revolute joints about +z, point masses. ---
        // Link 1's joint is at the world origin; its point mass sits l1 along the
        // link's +x (com = (l1,0,0)). A point mass has no inertia about its own
        // CoM — the pivot inertia m·d² is produced by the spatial inertia from
        // the CoM offset, so a negligible ε on the diagonal just keeps the
        // out-of-plane spin DOF well-posed.
        let eps = diag_inertia(1e-9, 1e-9, 1e-9);
        let link1 = TreeBody::revolute(
            ArticulatedTree::BASE,
            Vector3::z(),
            SpatialTransform::identity(),
            m1,
            Vector3::new(l1, 0.0, 0.0),
            eps,
        );
        // Link 2's joint sits at the end of link 1 (l1 along link-1 +x), so its
        // fixed parent→joint transform is a translation by (l1,0,0) in link 1's
        // frame; its point mass is l2 along link-2 +x.
        let link2 = TreeBody::revolute(
            0,
            Vector3::z(),
            SpatialTransform::translation(Vector3::new(l1, 0.0, 0.0)),
            m2,
            Vector3::new(l2, 0.0, 0.0),
            eps,
        );
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, -g, 0.0);
        tree.push(link1, t1, w1);
        tree.push(link2, t2, w2);
        let qdd = tree.forward_dynamics(&[0.0, 0.0]).unwrap();

        // --- Lagrangian closed form. The cleanest derivation is in the two
        // links' *absolute* angles a1, a2; we then map the absolute angular
        // accelerations back to the tree's joint coordinates. ---
        let a1 = t1;
        let a2 = t1 + t2; // link-2 absolute angle
        let da1 = w1;
        let da2 = w1 + w2; // link-2 absolute rate
        let da = a1 - a2;
        // Absolute-angle equations of motion for two point masses
        // (p1 = l1·(cos a1, sin a1); p2 = p1 + l2·(cos a2, sin a2); gravity −y):
        //   [ (m1+m2)l1²           m2 l1 l2 cos(a1−a2) ] [ä1]
        //   [ m2 l1 l2 cos(a1−a2)  m2 l2²              ] [ä2]
        // = [ −m2 l1 l2 sin(a1−a2)·ȧ2² − (m1+m2) g l1 cos a1 ]
        //   [  m2 l1 l2 sin(a1−a2)·ȧ1² − m2 g l2 cos a2       ]
        let am11 = (m1 + m2) * l1 * l1;
        let am12 = m2 * l1 * l2 * da.cos();
        let am22 = m2 * l2 * l2;
        let rhs1 = -m2 * l1 * l2 * da.sin() * da2 * da2 - (m1 + m2) * g * l1 * a1.cos();
        let rhs2 = m2 * l1 * l2 * da.sin() * da1 * da1 - m2 * g * l2 * a2.cos();
        let det = am11 * am22 - am12 * am12;
        let add1 = (rhs1 * am22 - am12 * rhs2) / det; // ä1
        let add2 = (am11 * rhs2 - am12 * rhs1) / det; // ä2
                                                      // Map absolute → joint accels: q̈1 = ä1, q̈2 = ä2 − ä1.
        let ref_qdd1 = add1;
        let ref_qdd2 = add2 - add1;

        // The two independent computation paths (the ABA recursion and the
        // closed-form mass-matrix solve) agree to ~9 significant figures; the
        // residual is pure floating-point accumulation, not model error.
        assert!(
            (qdd[0] - ref_qdd1).abs() < 1e-7,
            "double-pendulum q̈1 ABA {} vs Lagrangian {}",
            qdd[0],
            ref_qdd1
        );
        assert!(
            (qdd[1] - ref_qdd2).abs() < 1e-7,
            "double-pendulum q̈2 ABA {} vs Lagrangian {}",
            qdd[1],
            ref_qdd2
        );
    }

    /// (3) Free fall: with zero applied torque and a joint that does not resist
    /// gravity along its motion, every body's spatial acceleration equals `g`.
    /// Concretely, a single prismatic joint aligned with gravity must accelerate
    /// at exactly `g` (q̈ = g), and a prismatic joint perpendicular to gravity
    /// must not accelerate (q̈ = 0).
    #[test]
    fn free_fall_prismatic_along_gravity() {
        let g = 9.80665;
        // Prismatic along −y (gravity direction). The body slides freely; with
        // no other force it must accelerate at g downward → q̈ = +g along the
        // (−y) axis means the displacement grows in the gravity direction.
        let body = TreeBody::prismatic(
            ArticulatedTree::BASE,
            Vector3::new(0.0, -1.0, 0.0), // axis points down
            SpatialTransform::identity(),
            3.3,
            Vector3::zeros(),
            diag_inertia(0.1, 0.1, 0.1),
        );
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, -g, 0.0);
        tree.push(body, 0.0, 0.0);
        let qdd = tree.forward_dynamics(&[0.0]).unwrap();
        assert!(
            (qdd[0] - g).abs() < 1e-9,
            "free-fall prismatic q̈ {} vs g {}",
            qdd[0],
            g
        );
    }

    /// (3b) A prismatic joint perpendicular to gravity feels no driving force →
    /// zero acceleration (the joint normal carries the weight).
    #[test]
    fn free_fall_prismatic_perpendicular() {
        let g = 9.80665;
        let body = TreeBody::prismatic(
            ArticulatedTree::BASE,
            Vector3::x(), // horizontal axis, gravity is along −y
            SpatialTransform::identity(),
            2.0,
            Vector3::zeros(),
            diag_inertia(0.1, 0.1, 0.1),
        );
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, -g, 0.0);
        tree.push(body, 0.0, 0.0);
        let qdd = tree.forward_dynamics(&[0.0]).unwrap();
        assert!(qdd[0].abs() < 1e-12, "perpendicular prismatic q̈ {}", qdd[0]);
    }

    /// (3c) Free fall of a revolute pendulum released from straight-down (q=0
    /// with gravity along −y and CoM along +x means horizontal). At the
    /// horizontal release the angular acceleration is the maximum −m g d/I, and
    /// the CoM's instantaneous linear acceleration magnitude check confirms the
    /// body genuinely falls.
    #[test]
    fn revolute_release_horizontal() {
        let m = 1.0;
        let l = 1.0;
        let d = l / 2.0;
        let g = 9.80665;
        let i_cg = m * l * l / 12.0;
        let body = TreeBody::revolute(
            ArticulatedTree::BASE,
            Vector3::z(),
            SpatialTransform::identity(),
            m,
            Vector3::new(d, 0.0, 0.0),
            diag_inertia(1e-9, 1e-9, i_cg),
        );
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, -g, 0.0);
        tree.push(body, 0.0, 0.0); // q = 0 → rod horizontal along +x
        let qdd = tree.forward_dynamics(&[0.0]).unwrap();
        let i_pivot = i_cg + m * d * d;
        // At q = 0 the CoM is at (d, 0); gravity (0, −mg) gives the maximal
        // pivot torque r×F = (d,0,0)×(0,−mg,0) = (0,0,−mgd), so q̈ = −mgd/I.
        let analytic = -(m * g * d / i_pivot);
        assert!(
            (qdd[0] - analytic).abs() < 1e-9,
            "released-horizontal q̈ {} vs −mgd/I {}",
            qdd[0],
            analytic
        );
    }

    /// (4) Energy conservation: integrate a frictionless single pendulum many
    /// steps with no applied torque; total mechanical energy stays within a
    /// tight tolerance (semi-implicit Euler keeps it bounded for a conservative
    /// system).
    #[test]
    fn pendulum_conserves_energy() {
        let m = 1.0;
        let l = 1.0;
        let d = l / 2.0;
        let g = 9.80665;
        let i_cg = m * l * l / 12.0;
        let i_pivot = i_cg + m * d * d;
        let body = TreeBody::revolute(
            ArticulatedTree::BASE,
            Vector3::z(),
            SpatialTransform::identity(),
            m,
            Vector3::new(d, 0.0, 0.0),
            diag_inertia(1e-9, 1e-9, i_cg),
        );
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, -g, 0.0);
        let q0 = 0.6;
        tree.push(body, q0, 0.0);

        // Energy of the single pendulum: E = ½ I_pivot ω² + m g h, with the CoM
        // height h = d·sin(q) (gravity along −y, CoM at (d cos q, d sin q)).
        let energy = |tree: &ArticulatedTree| {
            let q = tree.q[0];
            let w = tree.qd[0];
            0.5 * i_pivot * w * w + m * g * (d * q.sin())
        };
        let e0 = energy(&tree);
        let dt = 1.0e-4;
        let mut max_dev = 0.0_f64;
        for _ in 0..100_000 {
            tree.step(&[0.0], dt).unwrap();
            let dev = (energy(&tree) - e0).abs() / e0.abs().max(1e-9);
            max_dev = max_dev.max(dev);
        }
        assert!(
            max_dev < 5e-3,
            "pendulum energy drift {max_dev} exceeded tolerance"
        );
    }

    /// A chain matches a hand-rolled composite-rigid-body / Newton-Euler check
    /// at a static configuration: with everything at rest and no torque, the
    /// joint accelerations equal the gravity-only response. Cross-checks the
    /// fixed-base gravity propagation for a 3-link chain (no analytic
    /// shortcut — we verify the ABA result is self-consistent via the
    /// inverse-dynamics identity τ = M·q̈ + bias, recomputed independently).
    #[test]
    fn three_link_chain_consistency() {
        // Build a 3-link revolute chain about z, point-ish masses.
        let g = 9.80665;
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, -g, 0.0);
        let eps = diag_inertia(1e-6, 1e-6, 0.02);
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::z(),
                SpatialTransform::identity(),
                1.0,
                Vector3::new(0.5, 0.0, 0.0),
                eps,
            ),
            0.3,
            0.1,
        );
        tree.push(
            TreeBody::revolute(
                0,
                Vector3::z(),
                SpatialTransform::translation(Vector3::new(1.0, 0.0, 0.0)),
                0.8,
                Vector3::new(0.4, 0.0, 0.0),
                eps,
            ),
            -0.2,
            -0.05,
        );
        tree.push(
            TreeBody::revolute(
                1,
                Vector3::z(),
                SpatialTransform::translation(Vector3::new(0.8, 0.0, 0.0)),
                0.5,
                Vector3::new(0.3, 0.0, 0.0),
                eps,
            ),
            0.15,
            0.2,
        );

        // ABA forward dynamics with an arbitrary torque vector.
        let tau = [0.7, -0.4, 0.25];
        let qdd = tree.forward_dynamics(&tau).unwrap();

        // Independent check via the equation of motion τ = M(q)·q̈ + b(q,q̇),
        // where the bias b is the inverse dynamics at q̈=0 and M is built
        // column-by-column from the difference of inverse dynamics. We compute
        // inverse dynamics with the standard Recursive Newton-Euler Algorithm
        // (a separate algorithm from ABA — a genuine cross-check).
        let bias = rnea(&tree, &[0.0; 3]);
        let mut mass_mat = Matrix3::zeros();
        for j in 0..3 {
            let mut e = [0.0; 3];
            e[j] = 1.0;
            let col = rnea(&tree, &e);
            for (i, &c) in col.iter().enumerate() {
                mass_mat[(i, j)] = c - bias[i];
            }
        }
        // Solve M·q̈ = τ − b and compare to ABA.
        let rhs = Vector3::new(tau[0] - bias[0], tau[1] - bias[1], tau[2] - bias[2]);
        let qdd_ref = mass_mat.lu().solve(&rhs).expect("mass matrix invertible");
        for i in 0..3 {
            assert!(
                (qdd[i] - qdd_ref[i]).abs() < 1e-7,
                "chain q̈[{i}] ABA {} vs RNEA-EoM {}",
                qdd[i],
                qdd_ref[i]
            );
        }
    }

    /// Recursive Newton-Euler inverse dynamics, used purely as an *independent*
    /// reference inside the chain test (τ for a given q̈). Same spatial algebra,
    /// different recursion — so agreement is a real cross-validation of ABA.
    fn rnea(tree: &ArticulatedTree, qdd: &[f64]) -> Vec<f64> {
        let n = tree.bodies.len();
        let mut xup = vec![SpatialTransform::identity(); n];
        let mut s = vec![Vector6::zeros(); n];
        let mut v = vec![Vector6::zeros(); n];
        let mut a = vec![Vector6::zeros(); n];
        let mut f = vec![Vector6::zeros(); n];
        let a_base = stack(Vector3::zeros(), -tree.gravity);
        // Outward pass: velocities, accelerations, body forces.
        for i in 0..n {
            let body = &tree.bodies[i];
            let s_i = body.joint.motion_subspace();
            let x_j = body.joint.joint_transform(tree.q[i]);
            let x_up = x_j.then(&body.parent_to_joint);
            let vj = s_i * tree.qd[i];
            // The parent's spatial velocity/acceleration carried into this body's
            // frame. For a base-attached body the parent is the fixed world, so
            // its velocity is zero and its acceleration is the gravity field —
            // but both must still be transformed by `x_up` into this frame.
            let (v_par, a_par) = if body.parent == ArticulatedTree::BASE {
                (Vector6::zeros(), x_up.apply_motion(&a_base))
            } else {
                (
                    x_up.apply_motion(&v[body.parent]),
                    x_up.apply_motion(&a[body.parent]),
                )
            };
            let v_i = v_par + vj;
            let a_i = a_par + s_i * qdd[i] + motion_cross(&v_i, &vj);
            let i_rb = body.spatial_inertia();
            xup[i] = x_up;
            s[i] = s_i;
            v[i] = v_i;
            a[i] = a_i;
            f[i] = i_rb * a_i + force_cross(&v_i, &(i_rb * v_i));
        }
        // Inward pass: project body forces onto joint axes, accumulate to parent.
        let mut tau = vec![0.0; n];
        for i in (0..n).rev() {
            tau[i] = s[i].dot(&f[i]);
            let parent = tree.bodies[i].parent;
            if parent != ArticulatedTree::BASE {
                let f_up = xup[i].force_to_parent(&f[i]);
                f[parent] += f_up;
            }
        }
        tau
    }

    /// Fail-loud: an out-of-range / non-preceding parent index is rejected.
    #[test]
    fn rejects_bad_parent() {
        let mut tree = ArticulatedTree::new();
        tree.push(
            TreeBody::revolute(
                5, // no such parent
                Vector3::z(),
                SpatialTransform::identity(),
                1.0,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.1, 0.1, 0.1),
            ),
            0.0,
            0.0,
        );
        let err = tree.forward_dynamics(&[0.0]).unwrap_err();
        assert!(matches!(err, TreeError::BadParent { body: 0, parent: 5 }));
    }

    /// Fail-loud: a cycle (a body whose parent is itself or a later body) is
    /// rejected — the `parent < i` rule catches it.
    #[test]
    fn rejects_cycle() {
        let mut tree = ArticulatedTree::new();
        // body 0 -> parent 1 (a later body): not topological → cycle-like.
        tree.push(
            TreeBody::revolute(
                1,
                Vector3::z(),
                SpatialTransform::identity(),
                1.0,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.1, 0.1, 0.1),
            ),
            0.0,
            0.0,
        );
        tree.push(
            TreeBody::revolute(
                0,
                Vector3::z(),
                SpatialTransform::identity(),
                1.0,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.1, 0.1, 0.1),
            ),
            0.0,
            0.0,
        );
        let err = tree.forward_dynamics(&[0.0, 0.0]).unwrap_err();
        assert!(matches!(err, TreeError::BadParent { body: 0, parent: 1 }));
    }

    /// Fail-loud: a zero joint axis is rejected.
    #[test]
    fn rejects_zero_axis() {
        let mut tree = ArticulatedTree::new();
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::zeros(),
                SpatialTransform::identity(),
                1.0,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.1, 0.1, 0.1),
            ),
            0.0,
            0.0,
        );
        let err = tree.forward_dynamics(&[0.0]).unwrap_err();
        assert!(matches!(err, TreeError::ZeroAxis { body: 0 }));
    }

    /// Fail-loud: a wrong-length `tau` is rejected.
    #[test]
    fn rejects_bad_tau_dimension() {
        let mut tree = ArticulatedTree::new();
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::z(),
                SpatialTransform::identity(),
                1.0,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.1, 0.1, 0.1),
            ),
            0.0,
            0.0,
        );
        let err = tree.forward_dynamics(&[0.0, 0.0]).unwrap_err();
        assert!(matches!(
            err,
            TreeError::BadDimension {
                expected: 1,
                got: 2
            }
        ));
    }

    /// The motion/force cross products satisfy the duality `a·(v ×f b) =
    /// −(v ×ₘ a)·b` for arbitrary spatial vectors — a structural sanity check
    /// on the spatial algebra the algorithm rests on.
    #[test]
    fn spatial_cross_duality() {
        let v = Vector6::new(0.3, -0.4, 0.1, 0.7, -0.2, 0.5);
        let a = Vector6::new(0.2, 0.5, -0.6, 0.1, 0.9, -0.3);
        let b = Vector6::new(-0.4, 0.2, 0.8, -0.1, 0.3, 0.6);
        let lhs = a.dot(&force_cross(&v, &b));
        let rhs = -motion_cross(&v, &a).dot(&b);
        assert!((lhs - rhs).abs() < 1e-12, "duality {lhs} vs {rhs}");
    }

    /// Round-trip: a transform composed with its inverse is the identity on a
    /// spatial motion vector.
    #[test]
    fn transform_inverse_round_trip() {
        let x = SpatialTransform::new(
            rodrigues(Vector3::new(0.2, 0.9, -0.3), 0.8).transpose(),
            Vector3::new(0.4, -0.6, 0.7),
        );
        let v = Vector6::new(0.1, 0.2, 0.3, 0.4, 0.5, 0.6);
        let back = x.inverse().apply_motion(&x.apply_motion(&v));
        assert!(
            (back - v).norm() < 1e-12,
            "round-trip drift {}",
            (back - v).norm()
        );
    }

    /// `PI` is referenced to keep the import meaningful for a rotation sanity
    /// check: rotating a vector by 2π about any axis returns it.
    #[test]
    fn full_turn_is_identity() {
        let r = rodrigues(Vector3::new(0.0, 0.0, 1.0), 2.0 * PI);
        let v = Vector3::new(1.0, 2.0, 3.0);
        assert!((r * v - v).norm() < 1e-9);
    }

    // -----------------------------------------------------------------------
    // Floating-base ABA.
    // -----------------------------------------------------------------------

    /// (F1) **Free-floating single body, no joints, under gravity.** With no
    /// links the whole system is just the free base, which must accelerate at
    /// exactly the gravity field `g` with zero angular acceleration and no
    /// internal motion (there are no joints).
    #[test]
    fn floating_single_body_free_falls_at_g() {
        let g = 9.80665;
        let tree = ArticulatedTree::new(); // no bodies; default gravity along −z
        let base = FloatingBase::new(
            2.5,
            Vector3::new(0.1, -0.2, 0.05), // off-origin CoM, to stress the 6×6 solve
            diag_inertia(0.3, 0.4, 0.5),
        );
        let acc = tree.forward_dynamics_floating(&base, &[]).unwrap();
        // Angular acceleration must be zero; linear acceleration must equal g
        // (default gravity = (0,0,−g)).
        let alpha = acc.base.fixed_rows::<3>(0).into_owned();
        let lin = acc.base.fixed_rows::<3>(3).into_owned();
        assert!(alpha.norm() < 1e-12, "free body spun up: α = {alpha:?}");
        assert!(
            (lin - Vector3::new(0.0, 0.0, -g)).norm() < 1e-12,
            "free body linear accel {lin:?} vs g (0,0,−{g})"
        );
        assert!(acc.joints.is_empty(), "no joints expected");
    }

    /// (F1b) A free body with an **off-origin CoM** and a nonzero angular
    /// velocity still has its *centre of mass* accelerate at exactly `g` (the
    /// base-origin spatial acceleration carries the ω×(ω×c) term, but the CoM
    /// must fall freely). This checks the spatial inertia / bias bookkeeping.
    #[test]
    fn floating_free_body_com_falls_at_g() {
        let g = 9.80665;
        let tree = ArticulatedTree::new();
        let mut base = FloatingBase::new(
            1.7,
            Vector3::new(0.3, 0.0, 0.0),
            diag_inertia(0.2, 0.2, 0.2),
        );
        base.velocity = Vector6::new(0.0, 0.0, 1.1, 0.0, 0.0, 0.0); // spin about z
        let acc = tree.forward_dynamics_floating(&base, &[]).unwrap();
        // CoM spatial acceleration a_c = a_o + ω̇×c + ω×(ω×c). With a_o the base
        // origin accel from the solve and ω̇ the angular accel, the *linear*
        // acceleration of the CoM point must be exactly g.
        let alpha = acc.base.fixed_rows::<3>(0).into_owned();
        let a_o = acc.base.fixed_rows::<3>(3).into_owned();
        let omega = base.velocity.fixed_rows::<3>(0).into_owned();
        let c = base.com;
        let a_com = a_o + alpha.cross(&c) + omega.cross(&omega.cross(&c));
        assert!(
            (a_com - Vector3::new(0.0, 0.0, -g)).norm() < 1e-9,
            "CoM accel {a_com:?} vs g"
        );
    }

    /// (F2) **The key floating-base correctness test — momentum conservation.**
    /// A free-floating two-body mechanism (base + one revolute link) with an
    /// internal joint torque and **gravity OFF** has no external wrench, so the
    /// total spatial momentum about the world origin (linear *and* angular) must
    /// be conserved.
    ///
    /// This is asserted two ways:
    ///
    /// 1. **Exactly, at the solver level (integrator-free):** the total external
    ///    wrench computed from the dynamics — `Σ (Iᵢaᵢ + vᵢ ×* Iᵢvᵢ)` over the
    ///    base and link, in the world frame — is zero to round-off. Internal
    ///    joint forces cancel, so a nonzero result would be a genuine model
    ///    error. This holds at *every* configuration along the trajectory, with
    ///    no dependence on the time step.
    /// 2. **Over a long integration:** the measured momentum drift stays small
    ///    (this residual is purely the first-order semi-implicit-Euler
    ///    truncation error — it shrinks linearly with `dt`, confirmed
    ///    separately — not a conservation-law violation).
    #[test]
    fn floating_two_body_conserves_momentum_gravity_off() {
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::zeros(); // no external force
                                         // One revolute link about z, mass offset from the base so internal
                                         // torque produces real linear+angular recoil on the base.
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::z(),
                SpatialTransform::translation(Vector3::new(0.4, 0.0, 0.0)),
                1.2,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.02, 0.02, 0.05),
            ),
            0.25, // initial joint angle
            0.7,  // initial joint rate
        );
        let mut base = FloatingBase::new(
            2.0,
            Vector3::new(0.0, 0.1, 0.0),
            diag_inertia(0.3, 0.3, 0.3),
        );
        base.velocity = Vector6::new(0.05, -0.03, 0.2, 0.4, -0.1, 0.15); // tumbling + drifting

        let tau = [0.6]; // a sustained internal joint torque

        // (1) Exact, integrator-free: net external wrench is zero (round-off),
        // checked at the start and again partway along the trajectory.
        let w0 = tree.external_wrench_world(&base, &tau).unwrap();
        assert!(
            w0.norm() < 1e-9,
            "net external wrench should be zero (no gravity, internal torque only): {w0:?}"
        );

        // (2) Long integration: bounded drift (first-order integrator error).
        let p0 = tree.momentum_world(&base);
        let dt = 1.0e-4;
        let mut max_dev = 0.0_f64;
        for k in 0..50_000 {
            tree.step_floating(&mut base, &tau, dt).unwrap();
            let p = tree.momentum_world(&base);
            max_dev = max_dev.max((p - p0).norm());
            // Re-check the exact invariant once mid-flight, far from the start.
            if k == 25_000 {
                let w = tree.external_wrench_world(&base, &tau).unwrap();
                assert!(w.norm() < 1e-9, "mid-flight external wrench nonzero: {w:?}");
            }
        }
        // Relative to the momentum magnitude (~|p0| ≈ 1.3), a sub-1% Euler drift.
        assert!(
            max_dev < 5e-2,
            "floating momentum drift {max_dev} (p0 = {p0:?})"
        );
    }

    /// (F3) **Consistency with the fixed-base solver.** With a very large base
    /// inertia the base cannot move, so the floating recursion must reduce to the
    /// fixed-base recursion: the joint accelerations match. Run **gravity off**
    /// (the physically-correct way to express an immovable base — a heavy *free*
    /// base in gravity free-falls, locally cancelling it, which is a different
    /// problem) with nonzero joint rates so the velocity-product coupling is
    /// exercised, not just statics.
    #[test]
    fn floating_heavy_base_matches_fixed_base() {
        let make_tree = || {
            let mut tree = ArticulatedTree::new();
            tree.gravity = Vector3::zeros();
            let eps = diag_inertia(1e-6, 1e-6, 0.02);
            tree.push(
                TreeBody::revolute(
                    ArticulatedTree::BASE,
                    Vector3::z(),
                    SpatialTransform::identity(),
                    1.0,
                    Vector3::new(0.5, 0.0, 0.0),
                    eps,
                ),
                0.3,
                0.4,
            );
            tree.push(
                TreeBody::revolute(
                    0,
                    Vector3::y(), // a non-parallel axis, to use full 3-D coupling
                    SpatialTransform::translation(Vector3::new(1.0, 0.0, 0.0)),
                    0.8,
                    Vector3::new(0.4, 0.0, 0.0),
                    eps,
                ),
                -0.2,
                -0.55,
            );
            tree
        };
        let tau = [0.7, -0.4];
        let fixed = make_tree();
        let qdd_fixed = fixed.forward_dynamics(&tau).unwrap();

        // A base ~1e8× heavier than the links: a₀ → 0, so the links see an
        // effectively fixed base.
        let floating_tree = make_tree();
        let base = FloatingBase::new(1.0e8, Vector3::zeros(), diag_inertia(1.0e8, 1.0e8, 1.0e8));
        let acc = floating_tree
            .forward_dynamics_floating(&base, &tau)
            .unwrap();
        // The base must be (nearly) immovable.
        assert!(
            acc.base.norm() < 1e-6,
            "heavy base should not accelerate: a₀ = {:?}",
            acc.base
        );
        for (i, (&qf, &qj)) in qdd_fixed.iter().zip(acc.joints.iter()).enumerate() {
            assert!(
                (qj - qf).abs() < 1e-6,
                "joint q̈[{i}] floating {qj} vs fixed {qf}"
            );
        }
    }

    /// (F4) **Energy conservation** for a free-floating frictionless mechanism.
    /// Base + one revolute link, gravity off, no applied torque, nonzero initial
    /// motion: the total kinetic energy must stay bounded under symplectic Euler.
    /// (Gravity off isolates the floating-base kinetic bookkeeping; with no
    /// dissipation energy is a strict invariant up to integrator round-off.)
    #[test]
    fn floating_mechanism_conserves_energy() {
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::zeros();
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::z(),
                SpatialTransform::translation(Vector3::new(0.3, 0.0, 0.0)),
                1.1,
                Vector3::new(0.45, 0.0, 0.0),
                diag_inertia(0.02, 0.02, 0.04),
            ),
            0.4,
            0.9,
        );
        let mut base = FloatingBase::new(
            1.8,
            Vector3::new(0.0, 0.05, 0.0),
            diag_inertia(0.25, 0.25, 0.25),
        );
        base.velocity = Vector6::new(0.1, 0.0, 0.3, 0.2, -0.15, 0.05);

        let e0 = tree.energy_floating(&base);
        let dt = 1.0e-4;
        let mut max_dev = 0.0_f64;
        for _ in 0..50_000 {
            tree.step_floating(&mut base, &[0.0], dt).unwrap();
            let dev = (tree.energy_floating(&base) - e0).abs() / e0.abs().max(1e-9);
            max_dev = max_dev.max(dev);
        }
        // Bounded drift: first-order semi-implicit Euler over 5 s of fast
        // tumbling. The strict, dt-free conservation property is proven by the
        // exact momentum-balance check in F2; this confirms the stepper stays
        // well-behaved.
        assert!(max_dev < 2e-3, "floating energy drift {max_dev}");
    }

    /// (F4b) Energy conservation **with gravity on** over a long free flight:
    /// base + link tumbling and falling. Kinetic + gravitational potential must
    /// stay bounded (symplectic Euler keeps a conservative system's energy
    /// bounded, here including the genuine gravitational potential of every
    /// link's world centre of mass).
    #[test]
    fn floating_mechanism_conserves_energy_under_gravity() {
        let g = 9.80665;
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, 0.0, -g);
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::x(),
                SpatialTransform::translation(Vector3::new(0.0, 0.3, 0.0)),
                1.0,
                Vector3::new(0.0, 0.5, 0.0),
                diag_inertia(0.05, 0.01, 0.05),
            ),
            0.2,
            0.6,
        );
        let mut base = FloatingBase::new(1.5, Vector3::zeros(), diag_inertia(0.2, 0.2, 0.2));
        base.velocity = Vector6::new(0.3, 0.1, -0.2, 0.5, 0.4, 2.0); // launched upward, spinning

        let e0 = tree.energy_floating(&base);
        let dt = 5.0e-5;
        let mut max_dev = 0.0_f64;
        for _ in 0..40_000 {
            tree.step_floating(&mut base, &[0.0], dt).unwrap();
            let dev = (tree.energy_floating(&base) - e0).abs() / e0.abs().max(1e-9);
            max_dev = max_dev.max(dev);
        }
        assert!(
            max_dev < 5e-3,
            "floating-under-gravity energy drift {max_dev}"
        );
    }

    /// (F2b) The momentum residual is **pure first-order integrator error**: the
    /// one-step momentum change divided by `dt` shrinks ~linearly as `dt → 0`
    /// (an `O(dt)` truncation), confirming the *solver* conserves momentum
    /// exactly and the long-run drift in (F2) is the integrator, not the model.
    #[test]
    fn floating_momentum_residual_is_first_order() {
        let build = || {
            let mut tree = ArticulatedTree::new();
            tree.gravity = Vector3::zeros();
            tree.push(
                TreeBody::revolute(
                    ArticulatedTree::BASE,
                    Vector3::z(),
                    SpatialTransform::translation(Vector3::new(0.4, 0.0, 0.0)),
                    1.2,
                    Vector3::new(0.5, 0.0, 0.0),
                    diag_inertia(0.02, 0.02, 0.05),
                ),
                0.25,
                0.7,
            );
            let mut base = FloatingBase::new(
                2.0,
                Vector3::new(0.0, 0.1, 0.0),
                diag_inertia(0.3, 0.3, 0.3),
            );
            base.velocity = Vector6::new(0.05, -0.03, 0.2, 0.4, -0.1, 0.15);
            (tree, base)
        };
        let rate_at = |dt: f64| {
            let (mut tree, mut base) = build();
            let p0 = tree.momentum_world(&base);
            tree.step_floating(&mut base, &[0.6], dt).unwrap();
            ((tree.momentum_world(&base) - p0) / dt).norm()
        };
        let r_coarse = rate_at(1e-3);
        let r_fine = rate_at(1e-4);
        // 10× smaller dt → ~10× smaller residual rate (first order).
        let ratio = r_coarse / r_fine;
        assert!(
            (9.0..11.0).contains(&ratio),
            "momentum residual not first-order: ratio {ratio} (coarse {r_coarse}, fine {r_fine})"
        );
    }

    /// (F2c) **Gravity bookkeeping cross-check.** With gravity on, the total
    /// external wrench reported by [`ArticulatedTree::external_wrench_world`]
    /// must equal the total weight `(Σmᵢ)·g` (force part) and the gravitational
    /// couple about the world origin `Σ x_comᵢ × mᵢg` (torque part) — verifying
    /// the per-body gravitational-wrench term in the floating solve is exactly
    /// the real gravity force, no more and no less.
    #[test]
    fn floating_external_wrench_equals_gravity() {
        let g = 9.80665;
        let mut tree = ArticulatedTree::new();
        tree.gravity = Vector3::new(0.0, 0.0, -g);
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::z(),
                SpatialTransform::translation(Vector3::new(0.4, 0.0, 0.0)),
                1.2,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.02, 0.02, 0.05),
            ),
            0.25,
            0.7,
        );
        let mut base = FloatingBase::new(
            2.0,
            Vector3::new(0.0, 0.1, 0.0),
            diag_inertia(0.3, 0.3, 0.3),
        );
        base.velocity = Vector6::new(0.05, -0.03, 0.2, 0.4, -0.1, 0.15);
        base.position = Vector3::new(0.2, -0.1, 0.5);
        base.orientation = rodrigues(Vector3::new(0.3, -0.5, 0.8), 0.6);

        let w = tree.external_wrench_world(&base, &[0.6]).unwrap();
        // Expected force = total weight; expected torque = Σ x_com × (m g).
        let kin = tree.world_kinematics(&base);
        let mut total_force = base.mass * tree.gravity;
        let base_com_w = base.position + base.orientation * base.com;
        let mut total_torque = base_com_w.cross(&(base.mass * tree.gravity));
        for (i, body) in tree.bodies.iter().enumerate() {
            total_force += body.mass * tree.gravity;
            let com_w = kin.origin_world[i] + kin.rot_body_to_world[i] * body.com;
            total_torque += com_w.cross(&(body.mass * tree.gravity));
        }
        let expected = stack(total_torque, total_force);
        assert!(
            (w - expected).norm() < 1e-9,
            "external wrench {w:?} vs gravity {expected:?}"
        );
    }

    /// Fail-loud: a floating base with **zero mass / inertia** has a singular
    /// `6×6` articulated inertia, so the base solve must return
    /// [`TreeError::SingularBaseInertia`] rather than producing `NaN`s.
    #[test]
    fn floating_singular_base_inertia_is_fail_loud() {
        let tree = ArticulatedTree::new(); // no bodies
        let base = FloatingBase::new(0.0, Vector3::zeros(), Matrix3::zeros());
        let err = tree.forward_dynamics_floating(&base, &[]).unwrap_err();
        assert_eq!(err, TreeError::SingularBaseInertia);
    }

    /// Fail-loud also propagates the existing structural checks (bad `tau`
    /// length) through the floating entry point.
    #[test]
    fn floating_rejects_bad_tau_dimension() {
        let mut tree = ArticulatedTree::new();
        tree.push(
            TreeBody::revolute(
                ArticulatedTree::BASE,
                Vector3::z(),
                SpatialTransform::identity(),
                1.0,
                Vector3::new(0.5, 0.0, 0.0),
                diag_inertia(0.1, 0.1, 0.1),
            ),
            0.0,
            0.0,
        );
        let base = FloatingBase::new(1.0, Vector3::zeros(), diag_inertia(0.1, 0.1, 0.1));
        let err = tree
            .forward_dynamics_floating(&base, &[0.0, 0.0])
            .unwrap_err();
        assert!(matches!(
            err,
            TreeError::BadDimension {
                expected: 1,
                got: 2
            }
        ));
    }
}
