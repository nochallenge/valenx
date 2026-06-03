//! Native in-process **transient structural dynamics** FEA solver —
//! Newmark-β implicit time integration.
//!
//! ## What this is
//!
//! A genuine, self-contained **transient (time-domain) dynamics**
//! solver for 4-node linear tetrahedra. The static solver
//! ([`crate::native_solver`]) answers "what is the deformed shape under
//! this load?"; the modal solver ([`crate::modal_solver`]) answers "at
//! what frequencies does the structure like to vibrate?". This module
//! answers the *time-history* question: **given an initial state and a
//! (possibly time-varying) load, how does the structure move, instant
//! by instant?**
//!
//! It solves the **semidiscrete equation of motion**
//!
//! ```text
//!   M·ü + C·u̇ + K·u = f(t)
//! ```
//!
//! — the mass `M` times acceleration, plus damping `C` times velocity,
//! plus stiffness `K` times displacement, balances the applied force.
//! `M` is the **consistent mass matrix** reused from
//! [`crate::modal_solver::assemble_global_stiffness_mass`]; `K` the same
//! elastic stiffness; `C` an optional **Rayleigh damping**
//! `C = α·M + β·K`.
//!
//! ## Newmark-β time integration
//!
//! The equation is marched in time by the **Newmark-β** method
//! (Newmark 1959) — the workhorse implicit integrator for structural
//! dynamics. Newmark assumes, over one step `Δt`,
//!
//! ```text
//!   u_{n+1} = u_n + Δt·u̇_n + Δt²·[(½−β)·ü_n + β·ü_{n+1}]
//!   u̇_{n+1} = u̇_n + Δt·[(1−γ)·ü_n + γ·ü_{n+1}]
//! ```
//!
//! Substituting these into the equation of motion at `t_{n+1}` gives a
//! linear system in the unknown `u_{n+1}` with the **effective
//! stiffness**
//!
//! ```text
//!   K_eff = K + (γ/(β·Δt))·C + (1/(β·Δt²))·M
//! ```
//!
//! and an effective load built from the previous step's state. Each
//! step is one solve of `K_eff·u_{n+1} = f_eff`; the velocity and
//! acceleration are then recovered from the Newmark relations. `K_eff`
//! is constant for a linear structure with a fixed step, so it is
//! factorised **once** and only back-substituted each step.
//!
//! With the standard **average-acceleration** parameters
//! `β = 1/4, γ = 1/2` ([`NewmarkParameters::average_acceleration`]) the
//! scheme is **unconditionally stable** and **second-order accurate**
//! and conserves energy for an undamped system — the time step is an
//! accuracy choice, not a stability bound.
//!
//! ## Honest scope
//!
//! This is a **real transient-dynamics v1** — the [tests](self) verify
//! that a single-DOF spring-mass oscillator reproduces its analytic
//! natural period `T = 2π·√(m/k)`. It is deliberately a v1:
//!
//! - **Linear** structural dynamics — `K` is the constant elastic
//!   stiffness. A geometrically- or materially-nonlinear transient
//!   solve would re-assemble `K` (and re-factorise `K_eff`) every step
//!   inside a Newton iteration.
//! - **Newmark-β** (average-acceleration by default). HHT-α (numerical
//!   damping of spurious high modes) is a one-parameter generalisation;
//!   an explicit central-difference scheme is the conditionally-stable
//!   alternative.
//! - **Rayleigh (mass + stiffness proportional) damping** only. Modal
//!   damping or a general `C` are follow-ups.
//! - **Tet4, isotropic, consistent mass** — inherited from
//!   [`crate::native_solver`] / [`crate::modal_solver`].
//!
//! Each of those is a documented, well-understood extension; none
//! affects the correctness of the linear transient response this
//! module computes.

use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};

use valenx_mesh::Mesh;

use crate::material::FemMaterial;
use crate::modal_solver::assemble_global_stiffness_mass;
use crate::native_solver::{ensure_finite3, NativeSolverError, NodalConstraint, NodalForce};

/// Round-8 hardening: cap on `controls.n_steps` accepted by
/// [`solve_transient_dynamics`]. Pre-fix, a hostile or runaway
/// caller passing `n_steps = 10_000_000_000` would either OOM the
/// host inside `Vec::with_capacity(controls.n_steps)` or silently
/// saturate-allocate. A million steps is well past any sensible
/// transient run (1M × default 1 ms Δt = 1000 s simulated time —
/// more than enough for the structural-dynamics use case).
pub const MAX_TIME_STEPS: usize = 1_000_000;

/// The two Newmark-β integration parameters `β` and `γ`.
///
/// Different `(β, γ)` choices give the classic Newmark family members.
/// [`NewmarkParameters::average_acceleration`] is the default and the
/// recommended choice for structural dynamics.
#[derive(Copy, Clone, Debug)]
pub struct NewmarkParameters {
    /// The Newmark `β` — controls how the acceleration is assumed to
    /// vary across the step. `1/4` (average acceleration) is
    /// unconditionally stable.
    pub beta: f64,
    /// The Newmark `γ` — controls the numerical damping. `1/2` gives
    /// no algorithmic damping (second-order accurate, energy
    /// conserving); `γ > 1/2` damps the high modes at the cost of
    /// dropping to first-order accuracy.
    pub gamma: f64,
}

impl NewmarkParameters {
    /// The **average-acceleration** (constant-average-acceleration,
    /// "trapezoidal") scheme — `β = 1/4, γ = 1/2`. Unconditionally
    /// stable, second-order accurate, energy-conserving for an
    /// undamped system. The standard choice.
    pub fn average_acceleration() -> NewmarkParameters {
        NewmarkParameters {
            beta: 0.25,
            gamma: 0.5,
        }
    }

    /// The **linear-acceleration** scheme — `β = 1/6, γ = 1/2`.
    /// Slightly more accurate than average-acceleration but only
    /// *conditionally* stable; offered for completeness.
    pub fn linear_acceleration() -> NewmarkParameters {
        NewmarkParameters {
            beta: 1.0 / 6.0,
            gamma: 0.5,
        }
    }
}

impl Default for NewmarkParameters {
    fn default() -> Self {
        NewmarkParameters::average_acceleration()
    }
}

/// An initial-condition entry for one node — a starting displacement
/// and/or velocity. A node with no [`NodalInitialState`] starts from
/// rest at the origin.
#[derive(Copy, Clone, Debug)]
pub struct NodalInitialState {
    /// 0-based node index.
    pub node: usize,
    /// Initial displacement `[ux, uy, uz]` (metres).
    pub displacement: [f64; 3],
    /// Initial velocity `[vx, vy, vz]` (m/s).
    pub velocity: [f64; 3],
}

/// Controls for a Newmark-β transient-dynamics run.
#[derive(Copy, Clone, Debug)]
pub struct DynamicsControls {
    /// The time step `Δt` (seconds). With the unconditionally-stable
    /// average-acceleration scheme this is an *accuracy* choice — a
    /// rule of thumb is `Δt ≤ T_min/20`, a twentieth of the period of
    /// the highest mode of interest.
    pub dt: f64,
    /// How many time steps to march.
    pub n_steps: usize,
    /// The Newmark parameters.
    pub newmark: NewmarkParameters,
    /// **Rayleigh mass-damping** coefficient `α` in `C = α·M + β·K`.
    /// `0` for an undamped run.
    pub rayleigh_alpha: f64,
    /// **Rayleigh stiffness-damping** coefficient `β` in
    /// `C = α·M + β·K`. `0` for an undamped run.
    pub rayleigh_beta: f64,
}

impl Default for DynamicsControls {
    /// Sensible defaults: `Δt = 1e-3 s`, 1000 steps, average-
    /// acceleration Newmark, undamped.
    fn default() -> Self {
        DynamicsControls {
            dt: 1.0e-3,
            n_steps: 1000,
            newmark: NewmarkParameters::average_acceleration(),
            rayleigh_alpha: 0.0,
            rayleigh_beta: 0.0,
        }
    }
}

/// One recorded instant of a transient-dynamics run.
#[derive(Clone, Debug)]
pub struct DynamicsFrame {
    /// Physical time `t` (seconds) of this frame.
    pub time: f64,
    /// Per-node displacement `[ux, uy, uz]` (metres) at this instant.
    pub displacement: Vec<[f64; 3]>,
}

/// Result of a transient-dynamics solve — the full time history.
#[derive(Clone, Debug)]
pub struct DynamicsSolution {
    /// One [`DynamicsFrame`] per time step, in time order. `frames[0]`
    /// is the state after the first step; the initial state is *not*
    /// stored as a frame.
    pub frames: Vec<DynamicsFrame>,
    /// The final physical time reached.
    pub final_time: f64,
}

impl DynamicsSolution {
    /// The displacement time history of one DOF (`node`, axis `0=x,
    /// 1=y, 2=z`) — `(time, displacement)` pairs, one per frame.
    ///
    /// Returns an empty vector if `node` / `axis` is out of range.
    pub fn dof_history(&self, node: usize, axis: usize) -> Vec<(f64, f64)> {
        if axis >= 3 {
            return Vec::new();
        }
        self.frames
            .iter()
            .filter_map(|f| {
                f.displacement.get(node).map(|d| (f.time, d[axis]))
            })
            .collect()
    }

    /// The peak absolute displacement reached by one DOF over the whole
    /// history — useful for reading a vibration amplitude.
    pub fn peak_abs_displacement(&self, node: usize, axis: usize) -> f64 {
        self.dof_history(node, axis)
            .iter()
            .map(|&(_, d)| d.abs())
            .fold(0.0, f64::max)
    }
}

/// Solve a linear transient-dynamics problem on a tetrahedral mesh by
/// Newmark-β implicit time integration.
///
/// `mesh` carries the Tet4 body; `material` supplies the elastic
/// constants **and the density** (dynamics needs mass). `constraints`
/// pin nodal DOFs (constrained DOFs are eliminated — held fixed for
/// all time). `forces` is the applied load, held **constant in time**
/// over the run (a step load applied at `t = 0`). `initial` sets the
/// starting displacement / velocity of any node (the rest start from
/// rest). `controls` sets the time step, the step count, the Newmark
/// parameters, and the Rayleigh damping.
///
/// Returns the [`DynamicsSolution`] — the displacement of every node
/// at every time step.
///
/// # Method
///
/// 1. Assemble the consistent mass `M` and stiffness `K`
///    ([`assemble_global_stiffness_mass`]); form the Rayleigh damping
///    `C = α·M + β·K`.
/// 2. Apply the constraints by **DOF elimination** (the correct
///    treatment for an integrator — a penalty spring would inject
///    spurious stiff modes).
/// 3. Recover the initial acceleration from `M·ü₀ = f − C·u̇₀ − K·u₀`.
/// 4. Form the constant Newmark effective stiffness `K_eff` and
///    factorise it once.
/// 5. Each step: build the effective load from the previous state,
///    back-substitute for `u_{n+1}`, and update the velocity and
///    acceleration from the Newmark relations.
///
/// # Errors
///
/// See [`NativeSolverError`] — empty mesh, no Tet4 block, a degenerate
/// element, bad connectivity, a bad material, or a factorisation
/// failure. A non-positive density yields
/// [`NativeSolverError::BadMaterial`].
pub fn solve_transient_dynamics(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    forces: &[NodalForce],
    initial: &[NodalInitialState],
    controls: &DynamicsControls,
) -> Result<DynamicsSolution, NativeSolverError> {
    // Round-8 guard: cap n_steps before any allocation. Pre-fix,
    // `Vec::with_capacity(controls.n_steps)` on a hostile n_steps
    // would either OOM the host or fall through to a silent
    // saturating allocation. MAX_TIME_STEPS = 1M is far past any
    // sensible transient run (a million 1ms steps = 1000 s simulated
    // time, more than enough for the structural-dynamics use case).
    if controls.n_steps > MAX_TIME_STEPS {
        return Err(NativeSolverError::InvalidParams {
            reason: format!(
                "n_steps {} exceeds the {MAX_TIME_STEPS}-step cap \
                 (DoS guard; lower the step count or coarsen Δt)",
                controls.n_steps
            ),
        });
    }
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    if !material.density.is_finite() || material.density <= 0.0 {
        return Err(NativeSolverError::BadMaterial(format!(
            "transient dynamics needs a positive density, got {}",
            material.density
        )));
    }
    // Round-1 H1–H3: dynamics had NO DOF cap — it walked straight into
    // the dense `assemble_global_stiffness_mass` (two `n_dof × n_dof`
    // matrices) and then built two more reduced dense matrices. Reject
    // an oversized mesh up front. (The assembler re-checks internally.)
    crate::native_solver::check_dense_dofs(n_nodes)?;

    // --- assemble the global K and consistent M ---
    let (k_full, m_full) = assemble_global_stiffness_mass(mesh, material)?;
    let n_dof = k_full.nrows();

    // --- which DOFs are free? constrained DOFs are eliminated ---
    let mut constrained = vec![false; n_dof];
    for c in constraints {
        if c.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: c.node,
                n_nodes,
            });
        }
        for (i, fixed) in c.fixed.iter().enumerate() {
            // A transient run constrains a DOF to a fixed value; only
            // homogeneous constraints (held at zero) are meaningful for
            // this v1, so any Some(_) eliminates the DOF.
            if fixed.is_some() {
                constrained[3 * c.node + i] = true;
            }
        }
    }
    let free: Vec<usize> = (0..n_dof).filter(|&d| !constrained[d]).collect();
    if free.is_empty() {
        // Nothing free to move — every frame is the rest state.
        // (Round-8: dropped a dead `.max(0)` on a usize here — usize
        // is always non-negative; the `.max(0)` was a leftover from
        // an earlier signed-int prototype.)
        let zero = vec![[0.0_f64; 3]; n_nodes];
        let frames = (0..controls.n_steps)
            .map(|s| DynamicsFrame {
                time: (s + 1) as f64 * controls.dt,
                displacement: zero.clone(),
            })
            .collect();
        return Ok(DynamicsSolution {
            frames,
            final_time: controls.n_steps as f64 * controls.dt,
        });
    }
    let n_free = free.len();

    // --- reduce K and M to the free DOFs (DOF elimination) ---
    let mut k = DMatrix::<f64>::zeros(n_free, n_free);
    let mut m = DMatrix::<f64>::zeros(n_free, n_free);
    for (ri, &gr) in free.iter().enumerate() {
        for (ci, &gc) in free.iter().enumerate() {
            k[(ri, ci)] = k_full[(gr, gc)];
            m[(ri, ci)] = m_full[(gr, gc)];
        }
    }
    // Rayleigh damping C = α·M + β·K.
    let c = &m * controls.rayleigh_alpha + &k * controls.rayleigh_beta;

    // --- the (time-constant) external load on the free DOFs ---
    let mut f_full = DVector::<f64>::zeros(n_dof);
    for force in forces {
        if force.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: force.node,
                n_nodes,
            });
        }
        // Round-1 H4: a non-finite force corrupts the effective load
        // every step — reject it before the run starts.
        ensure_finite3(&force.force, force.node, "force")?;
        for i in 0..3 {
            f_full[3 * force.node + i] += force.force[i];
        }
    }
    let mut f = DVector::<f64>::zeros(n_free);
    for (fi, &gd) in free.iter().enumerate() {
        f[fi] = f_full[gd];
    }

    // --- initial displacement / velocity on the free DOFs ---
    let mut u0_full = DVector::<f64>::zeros(n_dof);
    let mut v0_full = DVector::<f64>::zeros(n_dof);
    for init in initial {
        if init.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: init.node,
                n_nodes,
            });
        }
        // Round-1 H4: a non-finite initial displacement / velocity flows
        // into the initial-acceleration solve and then every Newmark
        // step, silently corrupting the whole history. Reject it.
        ensure_finite3(&init.displacement, init.node, "initial displacement")?;
        ensure_finite3(&init.velocity, init.node, "initial velocity")?;
        for i in 0..3 {
            u0_full[3 * init.node + i] = init.displacement[i];
            v0_full[3 * init.node + i] = init.velocity[i];
        }
    }
    let mut u = DVector::<f64>::zeros(n_free);
    let mut vel = DVector::<f64>::zeros(n_free);
    for (fi, &gd) in free.iter().enumerate() {
        u[fi] = u0_full[gd];
        vel[fi] = v0_full[gd];
    }

    // --- initial acceleration: M·a₀ = f − C·v₀ − K·u₀ ---
    let m_chol = m
        .clone()
        .cholesky()
        .ok_or(NativeSolverError::SolveFailed)?;
    let rhs0 = &f - &c * &vel - &k * &u;
    let mut acc = m_chol.solve(&rhs0);

    // --- the constant Newmark effective stiffness K_eff ---
    let dt = controls.dt.max(1e-30);
    let beta = controls.newmark.beta.max(1e-6);
    let gamma = controls.newmark.gamma;
    let a0 = 1.0 / (beta * dt * dt);
    let a1 = gamma / (beta * dt);
    let a2 = 1.0 / (beta * dt);
    let a3 = 1.0 / (2.0 * beta) - 1.0;
    let a4 = gamma / beta - 1.0;
    let a5 = dt * (gamma / (2.0 * beta) - 1.0);
    let a6 = dt * (1.0 - gamma);
    let a7 = dt * gamma;

    // K_eff = K + a0·M + a1·C.
    let k_eff = &k + &m * a0 + &c * a1;
    // Factorise K_eff once — it is constant for a linear structure with
    // a fixed step, so every step is just a back-substitution. K_eff is
    // SPD (K, M SPD; the C term keeps it so for the standard Newmark
    // parameters), so a sparse Cholesky is exact and fast. Convert the
    // dense K_eff to CSC.
    let k_eff_csc = dense_to_csc(&k_eff);
    let k_eff_chol = CscCholesky::factor(&k_eff_csc)
        .map_err(|_| NativeSolverError::SolveFailed)?;

    let mut frames = Vec::with_capacity(controls.n_steps);
    let mut time = 0.0;

    for _step in 0..controls.n_steps {
        // --- effective load at t_{n+1} ---
        // f_eff = f + M·(a0·u + a2·v + a3·a) + C·(a1·u + a4·v + a5·a).
        let m_term = &m * (&u * a0 + &vel * a2 + &acc * a3);
        let c_term = &c * (&u * a1 + &vel * a4 + &acc * a5);
        let f_eff = &f + m_term + c_term;

        // --- solve K_eff·u_{n+1} = f_eff ---
        let u_next_mat: DMatrix<f64> = k_eff_chol.solve(&f_eff);
        let u_next = u_next_mat.column(0).into_owned();

        // --- recover acceleration and velocity (Newmark relations) ---
        // a_{n+1} = a0·(u_{n+1}−u_n) − a2·v_n − a3·a_n.
        let acc_next =
            (&u_next - &u) * a0 - &vel * a2 - &acc * a3;
        // v_{n+1} = v_n + a6·a_n + a7·a_{n+1}.
        let vel_next = &vel + &acc * a6 + &acc_next * a7;

        u = u_next;
        vel = vel_next;
        acc = acc_next;
        time += dt;

        // --- record the frame (scatter the free DOFs back) ---
        let mut disp = vec![[0.0_f64; 3]; n_nodes];
        for (fi, &gd) in free.iter().enumerate() {
            disp[gd / 3][gd % 3] = u[fi];
        }
        // Constrained DOFs keep their (zero) prescribed value, already
        // the default in `disp`.
        frames.push(DynamicsFrame {
            time,
            displacement: disp,
        });
    }

    Ok(DynamicsSolution {
        final_time: time,
        frames,
    })
}

/// Convert a dense symmetric matrix to a CSC matrix (every non-zero
/// entry pushed once). Used to hand the Newmark effective stiffness to
/// the sparse Cholesky factoriser.
fn dense_to_csc(dense: &DMatrix<f64>) -> CscMatrix<f64> {
    let n = dense.nrows();
    let mut coo = CooMatrix::<f64>::new(n, dense.ncols());
    for j in 0..dense.ncols() {
        for i in 0..n {
            let v = dense[(i, j)];
            if v != 0.0 {
                coo.push(i, j, v);
            }
        }
    }
    CscMatrix::from(&coo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native_solver::structured_box_mesh;

    /// Node id of grid point `(i, j, k)` in a `structured_box_mesh`.
    fn nid(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
        i + (nx + 1) * j + (nx + 1) * (ny + 1) * k
    }

    #[test]
    fn newmark_average_acceleration_parameters() {
        // The default scheme is β=1/4, γ=1/2 — the unconditionally-
        // stable trapezoidal rule.
        let p = NewmarkParameters::average_acceleration();
        assert!((p.beta - 0.25).abs() < 1e-12);
        assert!((p.gamma - 0.5).abs() < 1e-12);
        assert_eq!(p.beta, NewmarkParameters::default().beta);
    }

    #[test]
    fn structure_at_rest_with_no_load_stays_at_rest() {
        // No initial motion, no force → every frame is the rest state.
        let mesh = structured_box_mesh(2.0, 1.0, 1.0, 4, 1, 1).expect("valid box params");
        let (nx, ny, nz) = (4, 1, 1);
        let mat = FemMaterial::default();
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let sol = solve_transient_dynamics(
            &mesh,
            &mat,
            &constraints,
            &[],
            &[],
            &DynamicsControls {
                dt: 1e-3,
                n_steps: 50,
                ..DynamicsControls::default()
            },
        )
        .unwrap();
        assert_eq!(sol.frames.len(), 50);
        for frame in &sol.frames {
            for d in &frame.displacement {
                let mag = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                assert!(mag < 1e-9, "a structure at rest must not move");
            }
        }
    }

    #[test]
    fn time_advances_by_dt_each_step() {
        // The physical clock advances by exactly Δt per step.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let mat = FemMaterial::default();
        let dt = 2.5e-4;
        let n = 20;
        let sol = solve_transient_dynamics(
            &mesh,
            &mat,
            &[NodalConstraint::fixed(0)],
            &[],
            &[],
            &DynamicsControls {
                dt,
                n_steps: n,
                ..DynamicsControls::default()
            },
        )
        .unwrap();
        assert_eq!(sol.frames.len(), n);
        assert!((sol.final_time - dt * n as f64).abs() < 1e-9);
        for (step, frame) in sol.frames.iter().enumerate() {
            assert!((frame.time - dt * (step + 1) as f64).abs() < 1e-9);
        }
    }

    #[test]
    fn spring_mass_oscillator_reproduces_its_natural_period() {
        // The headline transient-dynamics verification. A single-DOF
        // spring-mass oscillator vibrates freely with the analytic
        // period
        //   T = 2π·√(m/k) = 2π/ω.
        //
        // We make a finite-element model behave as exactly one
        // oscillator by **seeding the initial displacement with the
        // first mode shape itself**. A structure released from a pure
        // single-mode deflection (at rest) executes pure sinusoidal
        // motion in that one mode — `u(t) = φ₁·cos(ω₁·t)` — the exact
        // analogue of a single-DOF spring-mass oscillator. The Newmark
        // integration must then reproduce the period `T = 2π/ω₁`.
        //
        // The structure: a short bar along X, clamped at x=0. The modal
        // solver supplies both the first mode shape (the initial
        // condition) and the reference angular frequency ω₁.
        let (lx, ly, lz) = (1.0, 0.2, 0.2);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let mat = FemMaterial {
            youngs_modulus: 2.0e11,
            poisson_ratio: 0.3,
            density: 7800.0,
            ..FemMaterial::default()
        };

        // Clamp the x=0 face.
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }

        // --- the first mode shape + frequency from the modal solver ---
        // Seeding the transient run with mode 1 makes its free
        // vibration a single pure mode; ω₁ gives the analytic period.
        let modal = crate::modal_solver::solve_modal(
            &mesh,
            &mat,
            &constraints,
            1,
        )
        .unwrap();
        let omega = modal.modes[0].angular_frequency;
        assert!(omega > 0.0, "fundamental ω must be positive");
        let analytic_period = 2.0 * std::f64::consts::PI / omega;
        let mode1 = &modal.modes[0];

        // Initial displacement = the (scaled) first mode shape; released
        // from rest. The scale only sets the vibration amplitude.
        let amp_scale = 1.0e-3 / mode1.max_amplitude().max(1e-30);
        let initial: Vec<NodalInitialState> = (0..mesh.nodes.len())
            .map(|n| NodalInitialState {
                node: n,
                displacement: [
                    mode1.shape[n][0] * amp_scale,
                    mode1.shape[n][1] * amp_scale,
                    mode1.shape[n][2] * amp_scale,
                ],
                velocity: [0.0, 0.0, 0.0],
            })
            .collect();

        // March ≥ 4 full periods at ~80 steps per period.
        let dt = analytic_period / 80.0;
        let n_steps = (4.5 * 80.0) as usize;
        let sol = solve_transient_dynamics(
            &mesh,
            &mat,
            &constraints,
            &[],
            &initial,
            &DynamicsControls {
                dt,
                n_steps,
                newmark: NewmarkParameters::average_acceleration(),
                rayleigh_alpha: 0.0,
                rayleigh_beta: 0.0,
            },
        )
        .unwrap();
        assert_eq!(sol.frames.len(), n_steps);

        // Track the DOF that moves the most in the mode shape — its
        // history is a clean sinusoid `d₀·cos(ω·t)`.
        let mut track_node = 0;
        let mut track_axis = 0;
        let mut track_amp = 0.0_f64;
        for (n, s) in mode1.shape.iter().enumerate() {
            for (axis, &val) in s.iter().enumerate() {
                if val.abs() > track_amp {
                    track_amp = val.abs();
                    track_node = n;
                    track_axis = axis;
                }
            }
        }
        let d0 = track_amp * amp_scale; // initial value of the tracked DOF
        let history = sol.dof_history(track_node, track_axis);

        // The motion is a bounded, undamped oscillation about zero — it
        // neither grows nor decays (energy conserving Newmark).
        let peak = sol.peak_abs_displacement(track_node, track_axis);
        assert!(
            peak > 0.7 * d0 && peak < 1.3 * d0,
            "free vibration amplitude {peak} should stay near the initial {d0}"
        );

        // Measure the period from successive downward zero crossings of
        // the (single-mode → clean sinusoidal) history. Their spacing
        // is one full period.
        let mut crossings_down = Vec::new();
        for w in history.windows(2) {
            if w[0].1 > 0.0 && w[1].1 <= 0.0 {
                let frac = w[0].1 / (w[0].1 - w[1].1);
                crossings_down.push(w[0].0 + frac * (w[1].0 - w[0].0));
            }
        }
        assert!(
            crossings_down.len() >= 2,
            "the oscillator should complete multiple cycles, got {} crossings",
            crossings_down.len()
        );
        let measured_period = crossings_down[1] - crossings_down[0];
        let rel = (measured_period - analytic_period).abs() / analytic_period;
        assert!(
            rel < 0.03,
            "measured period {measured_period} should match analytic {analytic_period} (rel {rel})"
        );
        // Cross-check: every successive crossing pair is one period — a
        // single-mode oscillation is strictly periodic.
        for w in crossings_down.windows(2) {
            let t = w[1] - w[0];
            assert!(
                (t - analytic_period).abs() / analytic_period < 0.05,
                "every cycle should have the analytic period, got {t}"
            );
        }
    }

    #[test]
    fn rayleigh_damping_decays_the_vibration() {
        // With Rayleigh damping switched on, the free vibration of the
        // released bar must *decay* — the late-time amplitude is
        // smaller than the early-time amplitude. An undamped run does
        // not decay.
        let (lx, ly, lz) = (1.0, 0.2, 0.2);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let mat = FemMaterial {
            youngs_modulus: 2.0e11,
            poisson_ratio: 0.3,
            density: 7800.0,
            ..FemMaterial::default()
        };
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        // Release from the (scaled) first mode shape — a single-mode
        // initial condition, as in the natural-period test.
        let modal = crate::modal_solver::solve_modal(
            &mesh, &mat, &constraints, 1,
        )
        .unwrap();
        let mode1 = &modal.modes[0];
        let amp_scale = 1.0e-3 / mode1.max_amplitude().max(1e-30);
        let initial: Vec<NodalInitialState> = (0..mesh.nodes.len())
            .map(|n| NodalInitialState {
                node: n,
                displacement: [
                    mode1.shape[n][0] * amp_scale,
                    mode1.shape[n][1] * amp_scale,
                    mode1.shape[n][2] * amp_scale,
                ],
                velocity: [0.0, 0.0, 0.0],
            })
            .collect();
        // The DOF that moves most in the mode shape — tracked below.
        let mut track_node = 0;
        let mut track_axis = 0;
        let mut track_amp = 0.0_f64;
        for (n, s) in mode1.shape.iter().enumerate() {
            for (axis, &val) in s.iter().enumerate() {
                if val.abs() > track_amp {
                    track_amp = val.abs();
                    track_node = n;
                    track_axis = axis;
                }
            }
        }
        let period =
            2.0 * std::f64::consts::PI / modal.modes[0].angular_frequency;
        let dt = period / 80.0;
        let n_steps = (6.0 * 80.0) as usize;

        // Heavy mass-proportional damping.
        let sol = solve_transient_dynamics(
            &mesh,
            &mat,
            &constraints,
            &[],
            &initial,
            &DynamicsControls {
                dt,
                n_steps,
                newmark: NewmarkParameters::average_acceleration(),
                rayleigh_alpha: 80.0,
                rayleigh_beta: 0.0,
            },
        )
        .unwrap();
        let history = sol.dof_history(track_node, track_axis);
        // Amplitude in the first period vs the last period.
        let per_steps = 80;
        let early_peak = history
            .iter()
            .take(per_steps)
            .map(|&(_, d)| d.abs())
            .fold(0.0, f64::max);
        let late_peak = history
            .iter()
            .skip(history.len().saturating_sub(per_steps))
            .map(|&(_, d)| d.abs())
            .fold(0.0, f64::max);
        assert!(
            late_peak < 0.7 * early_peak,
            "damping should decay the vibration: early {early_peak}, late {late_peak}"
        );
    }

    #[test]
    fn rejects_bad_density() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let bad = FemMaterial {
            density: 0.0,
            ..FemMaterial::default()
        };
        let err = solve_transient_dynamics(
            &mesh,
            &bad,
            &[NodalConstraint::fixed(0)],
            &[],
            &[],
            &DynamicsControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::BadMaterial(_)));
    }

    #[test]
    fn rejects_a_mesh_without_tets() {
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = Mesh::new("surface");
        mesh.nodes.push(nalgebra::Vector3::new(0.0, 0.0, 0.0));
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 0, 0],
        });
        let err = solve_transient_dynamics(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[],
            &[],
            &DynamicsControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::NoTetBlock));
    }

    #[test]
    fn rejects_non_finite_initial_state() {
        // Round-1 H4: a NodalInitialState carrying a NaN displacement (or
        // velocity) flows into u0/v0 → the initial acceleration solve →
        // every Newmark step, silently corrupting the whole history.
        // Reject it up front.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let initial = vec![NodalInitialState {
            node: 7,
            displacement: [f64::NAN, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
        }];
        let err = solve_transient_dynamics(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[],
            &initial,
            &DynamicsControls {
                n_steps: 5,
                ..DynamicsControls::default()
            },
        )
        .unwrap_err();
        assert!(
            matches!(err, NativeSolverError::InvalidLoad { .. }),
            "expected InvalidLoad, got {err:?}"
        );

        // An infinite initial velocity is likewise rejected.
        let initial_v = vec![NodalInitialState {
            node: 7,
            displacement: [0.0, 0.0, 0.0],
            velocity: [0.0, f64::INFINITY, 0.0],
        }];
        let err_v = solve_transient_dynamics(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[],
            &initial_v,
            &DynamicsControls {
                n_steps: 5,
                ..DynamicsControls::default()
            },
        )
        .unwrap_err();
        assert!(matches!(err_v, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn rejects_nan_force() {
        // A NaN applied force corrupts the effective load each step.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let forces = vec![NodalForce {
            node: 7,
            force: [0.0, 0.0, f64::NAN],
        }];
        let err = solve_transient_dynamics(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &forces,
            &[],
            &DynamicsControls {
                n_steps: 5,
                ..DynamicsControls::default()
            },
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn rejects_n_steps_past_cap() {
        // Round-8 RED→GREEN: a runaway or hostile caller passing
        // `n_steps = 10_000_000_000` would either OOM the host inside
        // `Vec::with_capacity` or silently saturate-allocate. The
        // MAX_TIME_STEPS cap rejects the input up-front so the host
        // stays responsive.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let err = solve_transient_dynamics(
            &mesh,
            &FemMaterial::default(),
            &[NodalConstraint::fixed(0)],
            &[],
            &[],
            &DynamicsControls {
                n_steps: 10_000_000_000,
                ..DynamicsControls::default()
            },
        )
        .unwrap_err();
        match err {
            NativeSolverError::InvalidParams { reason } => {
                assert!(reason.contains("n_steps"), "msg: {reason}");
                assert!(
                    reason.contains(&MAX_TIME_STEPS.to_string()),
                    "msg: {reason}"
                );
            }
            other => panic!("expected InvalidParams, got {other:?}"),
        }
    }
}
