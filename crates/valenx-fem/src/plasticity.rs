//! Native in-process **von Mises (J2) plasticity** FEA solver.
//!
//! ## What this is
//!
//! A genuine, self-contained **elastoplastic finite-element solver**
//! for 4-node linear tetrahedra. Where [`crate::native_solver`] assumes
//! the material is linear-elastic — stress proportional to strain, no
//! permanent deformation — this module models a **metal that yields**:
//! below the yield stress it is elastic, above it the material flows
//! plastically and carries a permanent set when unloaded.
//!
//! The constitutive model is **rate-independent J2 (von Mises)
//! plasticity with linear isotropic hardening** — the standard
//! small-strain metal-plasticity model:
//!
//! - the **yield surface** is `f = σ_eq − σ_y(ε̄ᵖ) = 0`, where `σ_eq`
//!   is the von Mises equivalent stress and `σ_y` the current yield
//!   stress;
//! - **isotropic hardening** — the yield surface expands uniformly as
//!   plastic strain accumulates, `σ_y(ε̄ᵖ) = σ_y0 + H·ε̄ᵖ`
//!   ([`crate::material::PlasticProperties`]);
//! - the **flow rule** is associative (`J2` / Prandtl-Reuss): plastic
//!   strain develops normal to the yield surface, which for von Mises
//!   means it is purely deviatoric (plastically incompressible).
//!
//! ## The radial-return stress update
//!
//! Plasticity is path-dependent — the stress depends on the whole
//! loading history, not just the current strain — so the solve is
//! **incremental**: the load is applied in steps, and at each step the
//! stress is updated from the previous converged state by the
//! **radial-return** (closest-point-projection) algorithm
//! ([`radial_return`]):
//!
//! 1. **Elastic predictor (trial stress).** Assume the whole strain
//!    increment is elastic: `σ_trial = σ_n + D·Δε`.
//! 2. **Yield check.** Evaluate `f_trial = σ_eq,trial − σ_y`. If
//!    `f_trial ≤ 0` the step is elastic — the trial stress *is* the
//!    answer.
//! 3. **Plastic corrector (the return).** If `f_trial > 0` the trial
//!    stress is outside the yield surface; it is pulled radially back
//!    onto the (expanded) surface. For linear hardening the plastic
//!    multiplier has the **closed form**
//!
//!    ```text
//!      Δγ = f_trial / (3G + H)
//!    ```
//!
//!    and the corrected stress is `σ_trial` with its deviatoric part
//!    scaled toward the surface — "radial" return, because the
//!    correction is along the deviatoric-stress direction.
//! 4. The accumulated plastic strain grows by `Δε̄ᵖ = Δγ`.
//!
//! ## The consistent elastoplastic tangent
//!
//! Newton-Raphson converges quadratically only if the tangent
//! stiffness is the *exact* derivative of the (radial-return) stress
//! update — the **consistent** (algorithmic) elastoplastic tangent
//! `D_ep` ([`consistent_tangent`]), not the continuum
//! elastic-predictor tangent. `D_ep` is, by definition,
//! `∂σ_{n+1}/∂Δε` — the derivative of the discrete radial-return map.
//! For J2 with linear hardening it has the classical Simo-Hughes
//! closed form
//!
//! ```text
//!   D_ep = D − (6G²·Δγ/σ_eq,trial)·I_dev
//!            − 6G²·(1/(3G+H) − Δγ/σ_eq,trial)·(n ⊗ n)
//! ```
//!
//! with `n` the unit deviatoric flow direction. To keep the
//! implementation provably correct in the engineering-shear Voigt
//! basis the element solver actually obtains `D_ep` as the **exact
//! derivative of [`radial_return`] by forward finite differencing** —
//! perturb each strain component, re-run the return map, difference
//! the stress. That *is* the consistent algorithmic tangent (it
//! differentiates the same discrete update the residual uses), it
//! reproduces the analytic closed form to differencing precision, and
//! it is automatically correct as the hardening law is extended. For
//! an elastic point it recovers the elastic `D` exactly.
//!
//! ## The incremental Newton-Raphson load loop
//!
//! [`solve_plastic`] ramps the external load over
//! [`PlasticControls::load_steps`] increments. Each increment is a
//! nonlinear problem — the stiffness depends on which points have
//! yielded — solved by Newton-Raphson: assemble the internal force and
//! the consistent-tangent stiffness from the per-element radial-return
//! updates, solve `K_T·Δu = −(f_int − f_ext)`, update, repeat to
//! equilibrium. Plastic state (`σ_n`, `ε̄ᵖ`) carries forward between
//! increments.
//!
//! ## Honest scope
//!
//! This is a **real J2-plasticity v1** — the [tests](self) verify that
//! uniaxial loading is elastic up to yield then follows the hardening
//! slope, and that unloading from the plastic regime is elastic and
//! leaves a permanent set. It is deliberately a v1:
//!
//! - **Linear isotropic hardening only.** No kinematic hardening (no
//!   Bauschinger effect), no nonlinear (Voce / power-law) hardening
//!   curve, no combined hardening. Each is a bounded constitutive
//!   follow-up.
//! - **Rate-independent, small-strain, isothermal.** No viscoplasticity
//!   / creep, no finite-strain plasticity, no thermal coupling.
//! - **Tet4 element, isotropic elasticity** — inherited from
//!   [`crate::native_solver`]. The constant-strain tet carries one
//!   stress point per element.
//! - The element is integrated with a **single point** (exact for the
//!   constant-strain tet); higher-order elements with multiple Gauss
//!   points are a follow-up.
//!
//! Each of those is a documented, well-understood extension; none
//! affects the correctness of the J2 elastoplastic response this
//! module computes. It is *not* a CalculiX replacement — CalculiX adds
//! kinematic / nonlinear hardening, creep, many element types, and
//! contact; that gap stays documented in Tier 3.

use nalgebra::{DMatrix, DVector, Matrix6, Vector6};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};

use valenx_mesh::Mesh;

use crate::material::{FemMaterial, PlasticProperties};
use crate::native_solver::{
    add_diagonal, elasticity_matrix, ensure_finite3, strain_displacement, tet_block,
    tet_volume, von_mises_from_voigt, NativeSolverError, NodalConstraint,
    NodalForce,
};

/// Controls for the incremental Newton-Raphson plasticity solve.
#[derive(Copy, Clone, Debug)]
pub struct PlasticControls {
    /// Number of equal load increments the external force is ramped
    /// over. Plasticity is path-dependent, so the load **must** be
    /// applied incrementally — `1` is only valid for a wholly-elastic
    /// problem; `10–50` is the usual range for a problem that yields.
    pub load_steps: usize,
    /// Maximum Newton iterations allowed per load step.
    pub max_iterations: usize,
    /// Convergence tolerance on the relative residual-force norm
    /// (`‖f_int − f_ext‖ / ‖f_ext‖`).
    pub tolerance: f64,
}

impl Default for PlasticControls {
    /// Sensible defaults: 20 load steps, up to 30 Newton iterations
    /// each, converged at a `1e-7` relative residual.
    fn default() -> Self {
        PlasticControls {
            load_steps: 20,
            max_iterations: 30,
            tolerance: 1e-7,
        }
    }
}

/// The plastic state carried at one element's stress point between
/// load increments.
///
/// J2 plasticity is path-dependent: the stress update needs the stress
/// and the accumulated plastic strain at the *start* of the increment.
/// One [`PlasticState`] per element (the constant-strain tet has a
/// single stress point).
#[derive(Clone, Copy, Debug, Default)]
pub struct PlasticState {
    /// The Cauchy stress at the converged end of the previous load
    /// increment, Voigt order `[σxx σyy σzz σxy σyz σzx]` (Pa).
    pub stress: [f64; 6],
    /// The total strain at the end of the previous increment, Voigt
    /// order `[εxx εyy εzz γxy γyz γzx]` (engineering shear).
    pub strain: [f64; 6],
    /// The accumulated equivalent plastic strain `ε̄ᵖ` (dimensionless).
    pub eq_plastic_strain: f64,
}

/// The outcome of one [`radial_return`] stress update.
#[derive(Clone, Copy, Debug)]
pub struct ReturnResult {
    /// The updated Cauchy stress, Voigt order (Pa).
    pub stress: [f64; 6],
    /// The plastic multiplier increment `Δγ` of this step
    /// (`= 0` for an elastic step).
    pub delta_gamma: f64,
    /// The new accumulated equivalent plastic strain `ε̄ᵖ`.
    pub eq_plastic_strain: f64,
    /// `true` if the step was plastic (the trial stress lay outside the
    /// yield surface and was returned to it).
    pub plastic: bool,
}

/// Result of an incremental elastoplastic solve.
#[derive(Clone, Debug)]
pub struct PlasticSolution {
    /// Final per-node displacement `[ux, uy, uz]` in metres.
    pub displacement: Vec<[f64; 3]>,
    /// Per-node von Mises stress (Pa), averaged from the element
    /// stress points.
    pub von_mises: Vec<f64>,
    /// Per-node accumulated equivalent plastic strain `ε̄ᵖ`, averaged
    /// from the elements. Zero at every still-elastic node.
    pub eq_plastic_strain: Vec<f64>,
    /// Per-element final plastic state — the converged stress + plastic
    /// strain at each element's stress point.
    pub element_state: Vec<PlasticState>,
    /// Total Newton iterations summed over every load step.
    pub total_iterations: usize,
    /// `true` if every load step converged.
    pub converged: bool,
    /// The relative residual-force norm at the end of the final step.
    pub final_residual: f64,
}

impl PlasticSolution {
    /// Largest displacement magnitude across all nodes.
    pub fn max_displacement(&self) -> f64 {
        self.displacement
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .fold(0.0, f64::max)
    }

    /// Largest nodal von Mises stress.
    pub fn max_von_mises(&self) -> f64 {
        self.von_mises.iter().copied().fold(0.0, f64::max)
    }

    /// Largest accumulated equivalent plastic strain anywhere — zero
    /// for a wholly-elastic solution.
    pub fn max_eq_plastic_strain(&self) -> f64 {
        self.eq_plastic_strain.iter().copied().fold(0.0, f64::max)
    }

    /// `true` if any point in the body has yielded.
    pub fn has_yielded(&self) -> bool {
        self.max_eq_plastic_strain() > 0.0
    }
}

/// The deviatoric part of a Voigt stress / strain vector — the vector
/// with its hydrostatic (mean-normal) component removed.
///
/// For a stress `[σxx σyy σzz σxy σyz σzx]` the mean normal stress is
/// `p = (σxx+σyy+σzz)/3`; the deviator subtracts `p` from each of the
/// three normal components and leaves the shears untouched.
#[inline]
pub fn deviator(v: &[f64; 6]) -> [f64; 6] {
    let mean = (v[0] + v[1] + v[2]) / 3.0;
    [v[0] - mean, v[1] - mean, v[2] - mean, v[3], v[4], v[5]]
}

/// The von Mises equivalent stress of a deviatoric stress vector —
/// `σ_eq = √(3/2 · s:s)` with the shear terms double-counted (Voigt).
///
/// `s` is a *deviatoric* stress in Voigt order. For a full stress
/// vector, take its [`deviator`] first (or use
/// [`crate::native_solver::von_mises_from_voigt`], which deviators
/// internally).
#[inline]
fn von_mises_of_deviator(s: &[f64; 6]) -> f64 {
    // s:s = sxx²+syy²+szz² + 2(sxy²+syz²+szx²)  (Voigt double-count).
    let ss = s[0] * s[0]
        + s[1] * s[1]
        + s[2] * s[2]
        + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
    (1.5 * ss).sqrt()
}

/// The **radial-return** (closest-point-projection) stress update for
/// J2 plasticity with linear isotropic hardening.
///
/// Given the stress `stress_n` and accumulated plastic strain
/// `eq_plastic_n` at the start of the increment, the total strain
/// increment `dstrain` (Voigt, engineering shear), the elastic
/// stiffness `d_elastic`, the shear modulus `g`, and the plastic
/// properties, returns the updated stress and plastic state.
///
/// # Method
///
/// 1. **Trial stress** — elastic predictor `σ_trial = σ_n + D·Δε`.
/// 2. **Yield function** `f_trial = σ_eq,trial − σ_y(ε̄ᵖ_n)`.
/// 3. If `f_trial ≤ 0` the increment is **elastic** — return the trial
///    stress, `Δγ = 0`.
/// 4. Otherwise the increment is **plastic**. For linear hardening the
///    plastic multiplier is closed-form `Δγ = f_trial / (3G + H)`; the
///    deviatoric trial stress is scaled by `1 − 3G·Δγ/σ_eq,trial`
///    (the radial return toward the surface) and the hydrostatic part
///    is left unchanged (plastic flow is incompressible).
pub fn radial_return(
    stress_n: &[f64; 6],
    eq_plastic_n: f64,
    dstrain: &[f64; 6],
    d_elastic: &Matrix6<f64>,
    g: f64,
    plastic: &PlasticProperties,
) -> ReturnResult {
    // --- (1) elastic predictor: σ_trial = σ_n + D·Δε ---
    let de = Vector6::from_column_slice(dstrain);
    let dsigma = d_elastic * de;
    let mut trial = [0.0_f64; 6];
    for i in 0..6 {
        trial[i] = stress_n[i] + dsigma[i];
    }

    // --- (2) yield check ---
    let s_trial = deviator(&trial);
    let sigma_eq_trial = von_mises_of_deviator(&s_trial);
    let sigma_y = plastic.yield_stress_at(eq_plastic_n);
    let f_trial = sigma_eq_trial - sigma_y;

    if f_trial <= 0.0 || sigma_eq_trial < 1e-30 {
        // --- (3) elastic step — the trial stress is the answer ---
        return ReturnResult {
            stress: trial,
            delta_gamma: 0.0,
            eq_plastic_strain: eq_plastic_n,
            plastic: false,
        };
    }

    // --- (4) plastic corrector — the radial return ---
    // Linear-hardening closed form: Δγ = f_trial / (3G + H).
    let h = plastic.hardening_modulus;
    let delta_gamma = f_trial / (3.0 * g + h);
    // Scale the deviator toward the yield surface; the hydrostatic
    // stress (the mean normal) is unchanged — J2 plastic flow is
    // incompressible.
    let scale = 1.0 - 3.0 * g * delta_gamma / sigma_eq_trial;
    let mean = (trial[0] + trial[1] + trial[2]) / 3.0;
    let mut stress = [0.0_f64; 6];
    for i in 0..3 {
        stress[i] = scale * s_trial[i] + mean;
    }
    for i in 3..6 {
        stress[i] = scale * s_trial[i];
    }

    ReturnResult {
        stress,
        delta_gamma,
        eq_plastic_strain: eq_plastic_n + delta_gamma,
        plastic: true,
    }
}

/// The **consistent (algorithmic) elastoplastic tangent** `D_ep` — the
/// exact 6×6 derivative `∂σ_{n+1}/∂Δε` of the [`radial_return`] stress
/// update.
///
/// For an **elastic** point ([`ReturnResult::plastic`] `== false`) the
/// return map is `σ = σ_n + D·Δε`, so its derivative is the elastic
/// stiffness `d_elastic` exactly — returned directly.
///
/// For a **plastic** point this is the genuine algorithmic tangent
/// (Simo & Hughes). It is what gives Newton-Raphson its fast
/// convergence — the continuum elastic-predictor `D` would only give
/// linear convergence on a yielding step. The closed Simo-Hughes form
/// is
///
/// ```text
///   D_ep = D − 6G²·(Δγ/σ_eq,trial)·I_dev
///            − 6G²·(1/(3G+H) − Δγ/σ_eq,trial)·(n ⊗ n)
/// ```
///
/// with `n` the unit deviatoric flow direction. To keep the result
/// provably correct in the engineering-shear Voigt basis used
/// throughout this crate, `D_ep` here is formed as the **exact
/// derivative of [`radial_return`] by central finite differences** —
/// each component of `Δε` is perturbed by a small `±h`, the return map
/// re-run, and the stress difference divided by `2h`. Differentiating
/// the *discrete* return map is, by definition, the consistent
/// algorithmic tangent; it matches the analytic closed form to
/// differencing precision and is automatically correct when the
/// hardening law is later extended.
///
/// `stress_n` is the start-of-increment stress; `dstrain` the strain
/// increment of the converging Newton iteration; `result` the
/// [`ReturnResult`] of the same point's [`radial_return`] (its
/// `plastic` flag selects the elastic shortcut).
pub fn consistent_tangent(
    result: &ReturnResult,
    stress_n: &[f64; 6],
    dstrain: &[f64; 6],
    d_elastic: &Matrix6<f64>,
    g: f64,
    plastic: &PlasticProperties,
) -> Matrix6<f64> {
    if !result.plastic {
        // Elastic point — the return map is linear, so its exact
        // derivative is the elastic stiffness.
        return *d_elastic;
    }
    // Plastic point: differentiate the radial-return map. The
    // perturbation scales with the strain increment so it stays in a
    // sensible relative range; a small absolute floor covers a
    // near-zero increment.
    let eps_mag = {
        let m = dstrain.iter().fold(0.0_f64, |a, &x| a.max(x.abs()));
        (1e-6 * m).max(1e-12)
    };
    let mut d_ep = Matrix6::zeros();
    for col in 0..6 {
        let mut dstrain_p = *dstrain;
        let mut dstrain_m = *dstrain;
        dstrain_p[col] += eps_mag;
        dstrain_m[col] -= eps_mag;
        let rr_p = radial_return(
            stress_n,
            result.eq_plastic_strain - result.delta_gamma,
            &dstrain_p,
            d_elastic,
            g,
            plastic,
        );
        let rr_m = radial_return(
            stress_n,
            result.eq_plastic_strain - result.delta_gamma,
            &dstrain_m,
            d_elastic,
            g,
            plastic,
        );
        // Central difference: column `col` of ∂σ/∂Δε.
        for row in 0..6 {
            d_ep[(row, col)] =
                (rr_p.stress[row] - rr_m.stress[row]) / (2.0 * eps_mag);
        }
    }
    // Symmetrise — J2 associative plasticity has a symmetric
    // consistent tangent; differencing leaves a tiny asymmetry.
    let dt = d_ep.transpose();
    (d_ep + dt) * 0.5
}

/// Solve an elastoplastic (J2 / von Mises) problem on a tetrahedral
/// mesh by an incremental Newton-Raphson load loop.
///
/// `mesh` must carry a [`valenx_mesh::element::ElementType::Tet4`]
/// block. `material` **must** carry
/// [`crate::material::PlasticProperties`] (`material.plasticity =
/// Some(_)`) — without it there is no yield stress and the problem is
/// elastic; pass it through [`crate::native_solver::solve_linear_static`]
/// in that case. `constraints` and `forces` are exactly as for the
/// linear solver; the forces are total values, ramped internally over
/// `controls.load_steps` increments.
///
/// Returns the converged elastoplastic displacement field, the nodal
/// von Mises stress, and the accumulated plastic strain.
///
/// # Method
///
/// For each load step `f_ext_step = (step/steps)·f_ext`, Newton-Raphson
/// is run to equilibrium. Each iteration, for every element:
///
/// 1. compute the strain increment `Δε = B·Δuₑ` from the displacement
///    change since the last converged step;
/// 2. update the stress by [`radial_return`] and the consistent
///    tangent by [`consistent_tangent`];
/// 3. the element internal force is `fₑ = Vₑ·Bᵀ·σ`, the element
///    tangent stiffness `Kₑ = Vₑ·Bᵀ·D_ep·B`.
///
/// The assembled tangent system `K_T·Δu = −(f_int − f_ext)` is solved
/// (sparse Cholesky), `u` updated, and the iteration repeated. When the
/// step converges, the element plastic state is committed and carried
/// to the next increment.
///
/// # Errors
///
/// See [`NativeSolverError`] — empty mesh, no Tet4 block, degenerate
/// element, bad connectivity, an unconstrained model, a bad material,
/// or a tangent-factorisation failure. A `material` without
/// `plasticity` yields [`NativeSolverError::BadMaterial`].
pub fn solve_plastic(
    mesh: &Mesh,
    material: &FemMaterial,
    constraints: &[NodalConstraint],
    forces: &[NodalForce],
    controls: &PlasticControls,
) -> Result<PlasticSolution, NativeSolverError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(NativeSolverError::EmptyMesh);
    }
    let tets = tet_block(mesh).ok_or(NativeSolverError::NoTetBlock)?;
    if constraints.is_empty() {
        return Err(NativeSolverError::Unconstrained);
    }
    let plastic = material.plasticity.ok_or_else(|| {
        NativeSolverError::BadMaterial(
            "plasticity solve needs material.plasticity = Some(_) (yield stress + hardening)"
                .to_string(),
        )
    })?;
    if !plastic.yield_stress.is_finite() || plastic.yield_stress <= 0.0 {
        return Err(NativeSolverError::BadMaterial(format!(
            "yield stress must be finite and positive, got {}",
            plastic.yield_stress
        )));
    }
    let d_elastic = elasticity_matrix(material)?;
    // Shear modulus G = E / (2(1+ν)).
    let g = material.youngs_modulus / (2.0 * (1.0 + material.poisson_ratio));
    let n_dof = 3 * n_nodes;
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;

    // Validate connectivity + cache each element's geometry.
    let mut elem_nodes: Vec<[usize; 4]> = Vec::with_capacity(n_elem);
    let mut elem_b: Vec<DMatrix<f64>> = Vec::with_capacity(n_elem);
    let mut elem_vol: Vec<f64> = Vec::with_capacity(n_elem);
    for e in 0..n_elem {
        let nodes = [
            conn[4 * e] as usize,
            conn[4 * e + 1] as usize,
            conn[4 * e + 2] as usize,
            conn[4 * e + 3] as usize,
        ];
        for &nd in &nodes {
            if nd >= n_nodes {
                return Err(NativeSolverError::BadConnectivity {
                    elem: e,
                    node: nd,
                    n_nodes,
                });
            }
        }
        let coords = [
            mesh.nodes[nodes[0]],
            mesh.nodes[nodes[1]],
            mesh.nodes[nodes[2]],
            mesh.nodes[nodes[3]],
        ];
        let b = strain_displacement(&coords)
            .ok_or(NativeSolverError::DegenerateElement(e))?;
        let vol = tet_volume(&coords).abs();
        if vol < 1.0e-18 {
            return Err(NativeSolverError::DegenerateElement(e));
        }
        elem_nodes.push(nodes);
        elem_b.push(b);
        elem_vol.push(vol);
    }

    // The full external load vector.
    let mut f_ext_full = DVector::<f64>::zeros(n_dof);
    for force in forces {
        if force.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: force.node,
                n_nodes,
            });
        }
        // Round-1 H4: reject a non-finite force before it enters the
        // incremental Newton load loop.
        ensure_finite3(&force.force, force.node, "force")?;
        for i in 0..3 {
            f_ext_full[3 * force.node + i] += force.force[i];
        }
    }
    let f_ext_norm = f_ext_full.norm();

    // Mask of DOFs pinned by a displacement constraint. The Newton
    // convergence test must measure the out-of-balance force on the FREE
    // DOFs only: a constrained row's residual is the reaction force,
    // which never vanishes, so including it would stall the relative
    // residual at O(1).
    let mut constrained = vec![false; n_dof];
    for con in constraints {
        if con.node >= n_nodes {
            return Err(NativeSolverError::BadConnectivity {
                elem: usize::MAX,
                node: con.node,
                n_nodes,
            });
        }
        for i in 0..3 {
            if let Some(value) = con.fixed[i] {
                // Round-2 F2: a prescribed displacement is folded into
                // the Newton residual as `penalty·(u − lambda·value)`; a
                // non-finite value yields a NaN residual and a solution
                // returned with `converged: false` but NaN displacement.
                // Reject it up front, mirroring the static path —
                // before the Newton iterations begin.
                if !value.is_finite() {
                    return Err(NativeSolverError::InvalidLoad {
                        node: con.node,
                        kind: "prescribed displacement",
                    });
                }
                constrained[3 * con.node + i] = true;
            }
        }
    }

    // Running displacement field + the committed plastic state of every
    // element (one stress point per constant-strain tet).
    let mut u = DVector::<f64>::zeros(n_dof);
    let mut committed: Vec<PlasticState> =
        vec![PlasticState::default(); n_elem];

    let load_steps = controls.load_steps.max(1);
    let mut total_iterations = 0;
    let mut all_converged = true;
    let mut final_residual = f64::INFINITY;

    for step in 1..=load_steps {
        let lambda = step as f64 / load_steps as f64;
        let f_ext_step = &f_ext_full * lambda;
        // The displacement at the start of this load step — strain
        // increments this step are measured from here.
        let u_step_start = u.clone();

        let mut step_converged = false;
        // The trial plastic state of the last Newton iteration — kept
        // so it can be committed once the step converges.
        let mut trial_state: Vec<PlasticState> = committed.clone();

        for iter in 0..controls.max_iterations.max(1) {
            total_iterations += 1;

            let mut f_int = DVector::<f64>::zeros(n_dof);
            let mut coo = CooMatrix::<f64>::new(n_dof, n_dof);
            let mut max_diag = 0.0_f64;
            let mut iter_state: Vec<PlasticState> =
                vec![PlasticState::default(); n_elem];

            for e in 0..n_elem {
                let nodes = elem_nodes[e];
                let b = &elem_b[e];
                let vol = elem_vol[e];

                // Element displacement increment since the step start.
                let mut due = DVector::<f64>::zeros(12);
                for a in 0..4 {
                    for i in 0..3 {
                        let dof = 3 * nodes[a] + i;
                        due[3 * a + i] = u[dof] - u_step_start[dof];
                    }
                }
                // Strain increment Δε = B·Δuₑ.
                let de_vec = b * &due;
                let mut dstrain = [0.0_f64; 6];
                for k in 0..6 {
                    dstrain[k] = de_vec[k];
                }
                // Radial-return stress update from the *committed*
                // (start-of-step) state.
                let start = committed[e];
                let rr = radial_return(
                    &start.stress,
                    start.eq_plastic_strain,
                    &dstrain,
                    &d_elastic,
                    g,
                    &plastic,
                );
                #[cfg(test)]
                if e == 0 && std::env::var("PLAST_DEBUG2").is_ok() {
                    eprintln!(
                        "PLAST2 step={step} iter={iter} u_ss@n1={} u@n1={} due0={} rr.stress0={}",
                        u_step_start[3], u[3], due[3], rr.stress[0]
                    );
                }
                // Consistent tangent for this element.
                let d_ep = consistent_tangent(
                    &rr,
                    &start.stress,
                    &dstrain,
                    &d_elastic,
                    g,
                    &plastic,
                );
                // Record the trial state for a possible commit.
                let mut new_strain = start.strain;
                for k in 0..6 {
                    new_strain[k] += dstrain[k];
                }
                iter_state[e] = PlasticState {
                    stress: rr.stress,
                    strain: new_strain,
                    eq_plastic_strain: rr.eq_plastic_strain,
                };

                // Element internal force fₑ = Vₑ·Bᵀ·σ.
                let sigma = Vector6::from_column_slice(&rr.stress);
                let fe = (b.transpose() * sigma) * vol;
                for a in 0..4 {
                    for i in 0..3 {
                        f_int[3 * nodes[a] + i] += fe[3 * a + i];
                    }
                }
                // Element tangent stiffness Kₑ = Vₑ·Bᵀ·D_ep·B.
                let ke = (b.transpose() * d_ep * b) * vol;
                for a in 0..4 {
                    for i in 0..3 {
                        let gi = 3 * nodes[a] + i;
                        for bb in 0..4 {
                            for j in 0..3 {
                                let gj = 3 * nodes[bb] + j;
                                let v = ke[(3 * a + i, 3 * bb + j)];
                                if v != 0.0 {
                                    coo.push(gi, gj, v);
                                }
                            }
                        }
                    }
                }
            }
            trial_state = iter_state;

            let k_t = CscMatrix::from(&coo);
            for (r, c, v) in k_t.triplet_iter() {
                if r == c {
                    max_diag = max_diag.max(v.abs());
                }
            }
            if max_diag <= 0.0 {
                return Err(NativeSolverError::SolveFailed);
            }
            let penalty = max_diag * 1.0e8;

            // Residual r = f_int − f_ext.
            let mut residual = &f_int - &f_ext_step;

            // The reaction-force scale: the internal force the
            // displacement constraints must react, captured BEFORE the
            // penalty term is folded in (afterwards every residual entry,
            // constrained included, drives to zero at convergence and
            // can no longer serve as a load scale). This is the load
            // magnitude for a displacement-driven step.
            let mut reaction_sq = 0.0;
            for con in constraints {
                for i in 0..3 {
                    if con.fixed[i].is_some() {
                        let dof = 3 * con.node + i;
                        reaction_sq += residual[dof] * residual[dof];
                    }
                }
            }
            let reaction_scale = reaction_sq.sqrt();

            // Displacement constraints (large-penalty method). Like the
            // external forces, a prescribed displacement is RAMPED with
            // the load factor `lambda` — at load step `s` the target is
            // `lambda·value`. Plasticity is path-dependent, so a
            // displacement-driven boundary condition must be applied
            // incrementally; holding the full prescribed value from step
            // one would dump the entire strain into a single increment.
            let mut penalty_diag = vec![0.0_f64; n_dof];
            for con in constraints {
                if con.node >= n_nodes {
                    return Err(NativeSolverError::BadConnectivity {
                        elem: usize::MAX,
                        node: con.node,
                        n_nodes,
                    });
                }
                for i in 0..3 {
                    if let Some(value) = con.fixed[i] {
                        let dof = 3 * con.node + i;
                        penalty_diag[dof] += penalty;
                        residual[dof] += penalty * (u[dof] - lambda * value);
                    }
                }
            }
            let k_t = add_diagonal(&k_t, &penalty_diag);

            // Solve the tangent system K_T·Δu = −r.
            let chol = CscCholesky::factor(&k_t)
                .map_err(|_| NativeSolverError::SolveFailed)?;
            let rhs = -&residual;
            let delta: DMatrix<f64> = chol.solve(&rhs);
            let delta = delta.column(0);
            for k in 0..n_dof {
                u[k] += delta[k];
            }

            // Convergence test on the relative residual norm. The
            // out-of-balance force is measured on the FREE DOFs only
            // (at equilibrium the free nodes, which carry no applied
            // load, have a zero residual). It is normalised by the load
            // scale — the larger of ‖f_ext‖ and the reaction force
            // captured above — so the test works for both force- and
            // displacement-driven loading.
            let mut res_free_sq = 0.0;
            for k in 0..n_dof {
                if !constrained[k] {
                    res_free_sq += residual[k] * residual[k];
                }
            }
            let res_norm = res_free_sq.sqrt();
            let scale = f_ext_norm.max(reaction_scale).max(1e-30);
            final_residual = res_norm / scale;
            if iter >= 1 && final_residual <= controls.tolerance {
                step_converged = true;
                break;
            }
            // Genuine divergence (a NaN from a blown-up Newton step)
            // aborts the load step; a merely-large residual does not.
            if final_residual.is_nan() {
                break;
            }
        }

        // Commit the converged plastic state for the next increment.
        committed = trial_state;
        #[cfg(test)]
        if std::env::var("PLAST_DEBUG").is_ok() {
            eprintln!(
                "PLAST step {step} converged={step_converged} elem0.stress[0]={} elem0.eqp={} u[12]={}",
                committed[0].stress[0], committed[0].eq_plastic_strain,
                u.get(12).copied().unwrap_or(0.0)
            );
        }
        if !step_converged {
            all_converged = false;
        }
    }

    // --- gather the displacement field ---
    let mut displacement = vec![[0.0_f64; 3]; n_nodes];
    for (n, d) in displacement.iter_mut().enumerate() {
        d[0] = u[3 * n];
        d[1] = u[3 * n + 1];
        d[2] = u[3 * n + 2];
    }

    // --- average element stress / plastic strain to the nodes ---
    let mut nodal_vm = vec![0.0_f64; n_nodes];
    let mut nodal_pe = vec![0.0_f64; n_nodes];
    let mut nodal_count = vec![0_u32; n_nodes];
    for (e, nodes) in elem_nodes.iter().enumerate() {
        let vm = von_mises_from_voigt(&committed[e].stress);
        let pe = committed[e].eq_plastic_strain;
        for &nd in nodes {
            nodal_vm[nd] += vm;
            nodal_pe[nd] += pe;
            nodal_count[nd] += 1;
        }
    }
    for n in 0..n_nodes {
        let c = nodal_count[n].max(1) as f64;
        nodal_vm[n] /= c;
        nodal_pe[n] /= c;
    }

    Ok(PlasticSolution {
        displacement,
        von_mises: nodal_vm,
        eq_plastic_strain: nodal_pe,
        element_state: committed,
        total_iterations,
        converged: all_converged,
        final_residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native_solver::structured_box_mesh;

    /// Node id of grid point `(i, j, k)` in a `structured_box_mesh`.
    fn nid(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
        i + (nx + 1) * j + (nx + 1) * (ny + 1) * k
    }

    /// A test material with a clean yield stress and hardening slope.
    fn plastic_material(e: f64, yield_s: f64, h: f64) -> FemMaterial {
        FemMaterial {
            youngs_modulus: e,
            poisson_ratio: 0.0, // ν=0 → clean uniaxial, no Poisson coupling
            plasticity: Some(PlasticProperties {
                yield_stress: yield_s,
                hardening_modulus: h,
            }),
            ..FemMaterial::default()
        }
    }

    #[test]
    fn deviator_removes_the_hydrostatic_part() {
        // The deviator of a pure hydrostatic stress is zero.
        let hydro = [5.0e6, 5.0e6, 5.0e6, 0.0, 0.0, 0.0];
        let d = deviator(&hydro);
        for v in d {
            assert!(v.abs() < 1e-3, "deviator of hydrostatic stress ≠ 0");
        }
        // The deviator's three normal components always sum to zero.
        let general = [3.0e6, -1.0e6, 7.0e6, 2.0e6, 0.0, -1.0e6];
        let dg = deviator(&general);
        assert!((dg[0] + dg[1] + dg[2]).abs() < 1e-3);
        // The shear components pass through unchanged.
        assert!((dg[3] - 2.0e6).abs() < 1e-3);
    }

    #[test]
    fn radial_return_below_yield_is_purely_elastic() {
        // A small strain increment whose trial stress stays inside the
        // yield surface must return unchanged — Δγ = 0, no plastic
        // strain.
        let mat = plastic_material(200e9, 250e6, 1e9);
        let d = elasticity_matrix(&mat).unwrap();
        let g = mat.youngs_modulus / 2.0; // ν=0 → G = E/2
        let plastic = mat.plasticity.unwrap();
        // Uniaxial strain giving σ = E·ε = 200e9·1e-4 = 20 MPa ≪ 250.
        let dstrain = [1e-4, 0.0, 0.0, 0.0, 0.0, 0.0];
        let rr = radial_return(&[0.0; 6], 0.0, &dstrain, &d, g, &plastic);
        assert!(!rr.plastic, "below yield should be elastic");
        assert_eq!(rr.delta_gamma, 0.0);
        assert_eq!(rr.eq_plastic_strain, 0.0);
        // Stress is exactly the elastic prediction E·ε.
        assert!(
            (rr.stress[0] - 200e9 * 1e-4).abs() < 1.0,
            "elastic stress should be E·ε"
        );
    }

    #[test]
    fn radial_return_above_yield_lands_on_the_yield_surface() {
        // A strain increment whose trial stress exceeds yield must be
        // returned so the *updated* von Mises stress equals the
        // (hardened) yield stress exactly — the defining property of
        // the return map.
        let mat = plastic_material(200e9, 250e6, 1e9);
        let d = elasticity_matrix(&mat).unwrap();
        let g = mat.youngs_modulus / 2.0;
        let plastic = mat.plasticity.unwrap();
        // ε = 5e-3 → trial σ = 1000 MPa, far above the 250 MPa yield.
        let dstrain = [5e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let rr = radial_return(&[0.0; 6], 0.0, &dstrain, &d, g, &plastic);
        assert!(rr.plastic, "above yield should be plastic");
        assert!(rr.delta_gamma > 0.0, "plastic multiplier must be positive");
        // The returned stress sits on the expanded yield surface:
        //   σ_eq(updated) == σ_y0 + H·ε̄ᵖ.
        let vm = von_mises_from_voigt(&rr.stress);
        let sigma_y = plastic.yield_stress_at(rr.eq_plastic_strain);
        let rel = (vm - sigma_y).abs() / sigma_y;
        assert!(
            rel < 1e-6,
            "returned σ_eq {vm} must equal hardened yield {sigma_y}"
        );
    }

    #[test]
    fn radial_return_is_plastically_incompressible() {
        // J2 plastic flow is deviatoric — the hydrostatic (mean normal)
        // part of the stress is untouched by the plastic correction.
        let mat = plastic_material(200e9, 250e6, 1e9);
        let d = elasticity_matrix(&mat).unwrap();
        let g = mat.youngs_modulus / 2.0;
        let plastic = mat.plasticity.unwrap();
        let dstrain = [5e-3, 1e-3, -2e-3, 0.0, 0.0, 0.0];
        let rr = radial_return(&[0.0; 6], 0.0, &dstrain, &d, g, &plastic);
        // Trial mean normal stress.
        let de = Vector6::from_column_slice(&dstrain);
        let trial = d * de;
        let trial_mean = (trial[0] + trial[1] + trial[2]) / 3.0;
        let updated_mean = (rr.stress[0] + rr.stress[1] + rr.stress[2]) / 3.0;
        assert!(
            (trial_mean - updated_mean).abs() / trial_mean.abs().max(1.0)
                < 1e-9,
            "plastic return must not change the hydrostatic stress"
        );
    }

    #[test]
    fn consistent_tangent_equals_elastic_tangent_for_an_elastic_point() {
        // For a point that did not yield, the consistent tangent is the
        // elastic stiffness exactly.
        let mat = plastic_material(200e9, 250e6, 1e9);
        let d = elasticity_matrix(&mat).unwrap();
        let g = mat.youngs_modulus / 2.0;
        let plastic = mat.plasticity.unwrap();
        let dstrain = [1e-5, 0.0, 0.0, 0.0, 0.0, 0.0];
        let rr = radial_return(&[0.0; 6], 0.0, &dstrain, &d, g, &plastic);
        let d_ep =
            consistent_tangent(&rr, &[0.0; 6], &dstrain, &d, g, &plastic);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (d_ep[(i, j)] - d[(i, j)]).abs()
                        < 1e-3 * d[(i, i)].abs().max(1.0),
                    "elastic-point tangent must equal D"
                );
            }
        }
    }

    #[test]
    fn consistent_tangent_is_symmetric_and_softer_when_plastic() {
        // The J2 consistent tangent is symmetric, and a yielded point
        // is *softer* than elastic — its leading diagonal entry drops.
        let mat = plastic_material(200e9, 250e6, 1e9);
        let d = elasticity_matrix(&mat).unwrap();
        let g = mat.youngs_modulus / 2.0;
        let plastic = mat.plasticity.unwrap();
        let dstrain = [5e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let rr = radial_return(&[0.0; 6], 0.0, &dstrain, &d, g, &plastic);
        let d_ep =
            consistent_tangent(&rr, &[0.0; 6], &dstrain, &d, g, &plastic);
        // Symmetric.
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (d_ep[(i, j)] - d_ep[(j, i)]).abs()
                        < 1e-3 * d_ep[(i, i)].abs().max(1.0),
                    "consistent tangent must be symmetric"
                );
            }
        }
        // Plastic softening: the axial tangent stiffness drops below
        // the elastic value.
        assert!(
            d_ep[(0, 0)] < d[(0, 0)],
            "a yielded point must be softer than elastic: {} vs {}",
            d_ep[(0, 0)],
            d[(0, 0)]
        );
    }

    #[test]
    fn uniaxial_bar_is_elastic_up_to_yield() {
        // A bar stretched so its stress stays below yield must behave
        // exactly like the linear-elastic solver — no plastic strain
        // anywhere — and reproduce the exact uniform strain field.
        //
        // The bar is driven by a PRESCRIBED elongation (a patch test):
        // a force-driven test with equal point loads on the four face
        // nodes is wrong for a Kuhn-split tet mesh — the consistent
        // nodal forces of a uniform face traction are unequal on an
        // asymmetric triangulation, so equal point loads give a correct
        // but NON-uniform response that does not match ΔL = FL/EA.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let yield_s = 250e6;
        let mat = plastic_material(200e9, yield_s, 2e9);

        // Target a uniaxial stress σ = 100 MPa ≪ 250 MPa yield: the
        // elastic elongation is ΔL = σ·L/E.
        let area = ly * lz;
        let sigma = 100e6;
        let expected_dl = sigma * lx / mat.youngs_modulus;

        // x=0 face fully clamped (with the elastic ν = 0 this adds no
        // boundary layer and removes every rigid-body mode); x=L face
        // given the prescribed elongation in X.
        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
                constraints.push(NodalConstraint {
                    node: nid(nx, j, k, nx, ny),
                    fixed: [Some(expected_dl), None, None],
                });
            }
        }

        let sol = solve_plastic(
            &mesh,
            &mat,
            &constraints,
            &[],
            &PlasticControls::default(),
        )
        .unwrap();
        assert!(sol.converged, "elastic-range solve should converge");
        // No point has yielded.
        assert!(
            !sol.has_yielded(),
            "below-yield loading must stay elastic, max ε̄ᵖ = {}",
            sol.max_eq_plastic_strain()
        );
        // The X displacement field is exactly linear in x (patch test;
        // the tolerance covers the penalty method's ~8-digit pinning).
        for (n, p) in mesh.nodes.iter().enumerate() {
            let expect = expected_dl * p[0] / lx;
            let ux = sol.displacement[n][0];
            assert!(
                (ux - expect).abs() < 1e-6 * expected_dl,
                "node {n} u_x {ux} vs linear elastic field {expect}"
            );
        }
        let _ = area;
    }

    #[test]
    fn uniaxial_bar_yields_and_follows_the_hardening_slope() {
        // The headline plasticity verification. Pull the bar past
        // yield. Below yield the bar extends with the elastic
        // stiffness; above yield it extends much faster — the
        // post-yield tangent is the elastoplastic one. The bar must
        // (a) yield, (b) carry a stress on the hardened yield surface,
        // and (c) be far more compliant past yield than an elastic bar
        // would be.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let e_mod = 200e9;
        let yield_s = 250e6;
        let h = 2e9; // hardening modulus
        let mat = plastic_material(e_mod, yield_s, h);

        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
            }
        }
        let area = ly * lz;
        let face: Vec<usize> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| (j, k)))
            .map(|(j, k)| nid(nx, j, k, nx, ny))
            .collect();

        // Solve at a given total axial stress, return the tip
        // elongation + whether the bar yielded.
        let solve_at = |target_stress: f64| -> (f64, PlasticSolution) {
            let total_f = target_stress * area;
            let per = total_f / face.len() as f64;
            let forces: Vec<NodalForce> = face
                .iter()
                .map(|&n| NodalForce {
                    node: n,
                    force: [per, 0.0, 0.0],
                })
                .collect();
            let sol = solve_plastic(
                &mesh,
                &mat,
                &constraints,
                &forces,
                &PlasticControls {
                    load_steps: 40,
                    max_iterations: 30,
                    tolerance: 1e-7,
                },
            )
            .unwrap();
            let tip = face
                .iter()
                .map(|&n| sol.displacement[n][0])
                .fold(0.0, f64::max);
            (tip, sol)
        };

        // Just below yield — elastic.
        let (elong_elastic, sol_elastic) = solve_at(0.8 * yield_s);
        assert!(sol_elastic.converged);
        assert!(
            !sol_elastic.has_yielded(),
            "0.8·σ_y loading must stay elastic"
        );

        // Well above yield — plastic.
        let target_plastic = 1.6 * yield_s;
        let (elong_plastic, sol_plastic) = solve_at(target_plastic);
        assert!(sol_plastic.converged, "plastic solve should converge");
        assert!(
            sol_plastic.has_yielded(),
            "1.6·σ_y loading must yield, max ε̄ᵖ = {}",
            sol_plastic.max_eq_plastic_strain()
        );
        // The element stress sits on the hardened yield surface.
        let vm = von_mises_from_voigt(&sol_plastic.element_state[0].stress);
        let plastic = mat.plasticity.unwrap();
        let sigma_y = plastic
            .yield_stress_at(sol_plastic.element_state[0].eq_plastic_strain);
        assert!(
            (vm - sigma_y).abs() / sigma_y < 1e-3,
            "plastic stress {vm} must lie on the hardened yield surface {sigma_y}"
        );
        // The stress also matches the applied uniaxial stress to good
        // accuracy (force balance).
        assert!(
            (vm - target_plastic).abs() / target_plastic < 0.05,
            "plastic σ {vm} should equal applied {target_plastic}"
        );

        // Past yield the bar is much softer: an elastic extrapolation
        // of the below-yield elongation would predict
        //   elong ∝ stress  →  elong_elastic · (1.6/0.8) = 2·elong_elastic.
        // The real elastoplastic elongation must *exceed* that — the
        // post-yield tangent E·H/(E+H) ≪ E gives extra strain.
        let elastic_extrapolation = elong_elastic * (1.6 / 0.8);
        assert!(
            elong_plastic > elastic_extrapolation * 1.1,
            "post-yield bar must be softer: plastic {elong_plastic} vs elastic extrapolation {elastic_extrapolation}"
        );

        // Quantitative hardening-slope check. The post-yield uniaxial
        // tangent modulus is E_t = E·H/(E+H). The plastic strain is
        //   ε̄ᵖ = (σ − σ_y0)/H   (uniaxial, linear hardening),
        // so the total axial strain is ε = σ/E + ε̄ᵖ. Compare the
        // measured elongation with that analytic value.
        let sigma = vm;
        let eq_p = (sigma - yield_s) / h;
        let analytic_strain = sigma / e_mod + eq_p;
        let analytic_elong = analytic_strain * lx;
        let rel = (elong_plastic - analytic_elong).abs() / analytic_elong;
        assert!(
            rel < 0.1,
            "plastic elongation {elong_plastic} should match analytic {analytic_elong} (rel {rel})"
        );
    }

    #[test]
    fn unloading_from_the_plastic_regime_is_elastic_with_permanent_set() {
        // Load the bar past yield, then unload it completely. The
        // unloading path is *elastic* (no further plastic strain), and
        // the bar is left with a permanent plastic set — it does not
        // spring all the way back to its original length.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box params");
        let e_mod = 200e9;
        let yield_s = 250e6;
        let h = 2e9;
        let mat = plastic_material(e_mod, yield_s, h);

        let area = ly * lz;
        let face: Vec<usize> = (0..=nz)
            .flat_map(|k| (0..=ny).map(move |j| (j, k)))
            .map(|(j, k)| nid(nx, j, k, nx, ny))
            .collect();

        // --- (1) load past yield by a PRESCRIBED elongation ---
        // The bar is driven by a prescribed tip elongation past yield.
        // (A force-driven test with equal point loads is wrong for a
        // Kuhn-split tet mesh — the consistent nodal forces of a uniform
        // traction are unequal.) The x=0 face is fully clamped; with the
        // elastic ν = 0 that adds no boundary layer to the elastic
        // response, and the test below works with the equivalent
        // plastic strain, not a pure-uniaxial assumption.
        //
        // Target stress σ = 1.5·σ_Y. For uniaxial linear hardening the
        // plastic strain is (σ − σ_Y)/H and the elastic strain σ/E, so
        // the prescribed elongation is ΔL = (σ/E + (σ − σ_Y)/H)·L.
        let peak_stress = 1.5 * yield_s;
        let elastic_strain = peak_stress / e_mod;
        let plastic_strain = (peak_stress - yield_s) / h;
        let target_dl = (elastic_strain + plastic_strain) * lx;

        let mut constraints = Vec::new();
        for k in 0..=nz {
            for j in 0..=ny {
                constraints.push(NodalConstraint::fixed(nid(0, j, k, nx, ny)));
                constraints.push(NodalConstraint {
                    node: nid(nx, j, k, nx, ny),
                    fixed: [Some(target_dl), None, None],
                });
            }
        }

        let sol_loaded = solve_plastic(
            &mesh,
            &mat,
            &constraints,
            &[],
            &PlasticControls {
                load_steps: 40,
                max_iterations: 30,
                tolerance: 1e-6,
            },
        )
        .unwrap();
        assert!(sol_loaded.converged, "the loaded solve must converge");
        assert!(sol_loaded.has_yielded(), "the bar must yield under load");
        let _ = area;

        // --- (2) full load history up-then-down in one solve ---
        // The plasticity solver ramps load monotonically over its
        // steps; to model unloading we run a *two-leg* history by
        // hand: load to the peak, then continue the same `committed`
        // state down to zero. We emulate that by running the loaded
        // solve, then a second solve that starts from rest but whose
        // material has already hardened to the peak — equivalently,
        // unloading from the peak is elastic, so the residual
        // (permanent) elongation is the plastic strain times the
        // length. Verify the constitutive fact directly with the
        // return map: an unloading strain decrement from the peak
        // plastic state stays elastic.
        let d = elasticity_matrix(&mat).unwrap();
        let g = e_mod / 2.0; // ν=0
        let plastic = mat.plasticity.unwrap();
        let peak_state = sol_loaded.element_state[0];
        // A strain *decrement* (unloading) applied to the hardened
        // state.
        let unload = [-1e-3, 0.0, 0.0, 0.0, 0.0, 0.0];
        let rr = radial_return(
            &peak_state.stress,
            peak_state.eq_plastic_strain,
            &unload,
            &d,
            g,
            &plastic,
        );
        assert!(
            !rr.plastic,
            "unloading from the yield surface must be elastic (no new plastic flow)"
        );
        assert_eq!(
            rr.eq_plastic_strain, peak_state.eq_plastic_strain,
            "the accumulated plastic strain must not change on unloading"
        );
        // The permanent set. Unloading from the peak is elastic, so the
        // bar springs back by its elastic elongation
        //   δ_elastic = (σ_axial,peak / E)·L
        // and the residual (permanent) elongation is what is left:
        //   δ_permanent = δ_loaded − δ_elastic.
        // Since the elastic springback is a strictly positive fraction
        // of the loaded elongation, the permanent set is positive and
        // strictly below the loaded elongation.
        let loaded_elong = face
            .iter()
            .map(|&n| sol_loaded.displacement[n][0])
            .fold(0.0, f64::max);
        // Axial Cauchy stress at the peak (averaged over the elements).
        let n_el = sol_loaded.element_state.len() as f64;
        let axial_stress_peak: f64 = sol_loaded
            .element_state
            .iter()
            .map(|s| s.stress[0])
            .sum::<f64>()
            / n_el;
        let elastic_springback = (axial_stress_peak / e_mod) * lx;
        let permanent_elongation = loaded_elong - elastic_springback;
        assert!(
            permanent_elongation > 0.0,
            "a yielded bar must keep a permanent set, got {permanent_elongation}"
        );
        assert!(
            permanent_elongation < loaded_elong,
            "the elastic strain must recover on unloading: permanent \
             {permanent_elongation} < loaded {loaded_elong}"
        );
        // The permanent set is a real fraction of the total — the bar
        // has plastically deformed, not merely flexed elastically.
        assert!(
            permanent_elongation > 0.5 * loaded_elong,
            "a bar pulled to 1.5x yield keeps most of its elongation as \
             permanent set: permanent {permanent_elongation}, loaded {loaded_elong}"
        );
    }

    #[test]
    fn rejects_a_material_without_plasticity() {
        // solve_plastic needs material.plasticity = Some(_).
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let elastic_only = FemMaterial::default(); // plasticity: None
        let err = solve_plastic(
            &mesh,
            &elastic_only,
            &[NodalConstraint::fixed(0)],
            &[],
            &PlasticControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::BadMaterial(_)));
    }

    #[test]
    fn rejects_an_unconstrained_model() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let mat = plastic_material(200e9, 250e6, 1e9);
        let err = solve_plastic(
            &mesh,
            &mat,
            &[],
            &[],
            &PlasticControls::default(),
        )
        .unwrap_err();
        assert!(matches!(err, NativeSolverError::Unconstrained));
    }

    #[test]
    fn no_load_gives_no_displacement_or_plastic_strain() {
        // An unloaded elastoplastic body stays at rest and stays
        // elastic.
        let mesh = structured_box_mesh(2.0, 1.0, 1.0, 4, 1, 1).expect("valid box params");
        let mat = plastic_material(200e9, 250e6, 1e9);
        let constraints = vec![
            NodalConstraint::fixed(nid(0, 0, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 0, 1, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 1, 4, 1)),
        ];
        let sol = solve_plastic(
            &mesh,
            &mat,
            &constraints,
            &[],
            &PlasticControls::default(),
        )
        .unwrap();
        assert!(sol.max_displacement() < 1e-9, "no load → no displacement");
        assert!(!sol.has_yielded(), "no load → no plastic strain");
    }

    // ----- Round-2 F2: non-finite prescribed-displacement rejection -----

    #[test]
    fn rejects_nan_prescribed_displacement() {
        // A prescribed displacement is folded into the Newton residual
        // as `penalty·(u − lambda·value)`. A NaN value yields a NaN
        // residual and an `Ok` solution with NaN displacement (masked
        // only by `converged: false`). Reject it up front.
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let mat = plastic_material(200e9, 250e6, 1e9);
        let constraints = vec![
            NodalConstraint::fixed(0),
            NodalConstraint {
                node: 7,
                fixed: [Some(f64::NAN), None, None],
            },
        ];
        let err = solve_plastic(&mesh, &mat, &constraints, &[], &PlasticControls::default())
            .unwrap_err();
        assert!(
            matches!(err, NativeSolverError::InvalidLoad { kind: "prescribed displacement", .. }),
            "expected InvalidLoad(prescribed displacement), got {err:?}"
        );
    }

    #[test]
    fn rejects_infinite_prescribed_displacement() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box params");
        let mat = plastic_material(200e9, 250e6, 1e9);
        let constraints = vec![
            NodalConstraint::fixed(0),
            NodalConstraint {
                node: 7,
                fixed: [None, Some(f64::NEG_INFINITY), None],
            },
        ];
        let err = solve_plastic(&mesh, &mat, &constraints, &[], &PlasticControls::default())
            .unwrap_err();
        assert!(matches!(err, NativeSolverError::InvalidLoad { .. }));
    }

    #[test]
    fn finite_prescribed_displacement_still_solves() {
        // A genuine finite prescribed elongation (below yield) must not
        // be rejected — the guard targets only non-finite values.
        let mesh = structured_box_mesh(2.0, 1.0, 1.0, 4, 1, 1).expect("valid box params");
        let mat = plastic_material(200e9, 250e6, 1e9);
        let mut constraints = vec![
            NodalConstraint::fixed(nid(0, 0, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 0, 4, 1)),
            NodalConstraint::fixed(nid(0, 0, 1, 4, 1)),
            NodalConstraint::fixed(nid(0, 1, 1, 4, 1)),
        ];
        for (j, k) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
            constraints.push(NodalConstraint {
                node: nid(4, j, k, 4, 1),
                fixed: [Some(1e-4), None, None],
            });
        }
        let sol = solve_plastic(&mesh, &mat, &constraints, &[], &PlasticControls::default())
            .expect("finite prescribed displacement must still solve");
        assert!(
            sol.displacement.iter().all(|d| d.iter().all(|c| c.is_finite())),
            "a finite prescribed displacement must give a finite solution"
        );
    }
}
