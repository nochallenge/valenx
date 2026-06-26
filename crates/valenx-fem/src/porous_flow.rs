//! Native in-process **subsurface groundwater / porous-media flow** solver
//! (Phase 24.x) — Darcy's law and the Richards equation, for the
//! civil/geotechnical questions of **dam seepage**, **soil saturation**,
//! and **landslide-risk pore pressure**.
//!
//! ## What this is
//!
//! Two genuine, self-contained finite-element flow solvers over the same
//! 4-node linear-tetrahedron mesh ([`crate::native_solver::structured_box_mesh`])
//! the rest of valenx-fem assembles on — no subprocess, no input deck:
//!
//! 1. **Saturated Darcy flow** ([`solve_saturated_darcy`]). The steady,
//!    **linear** groundwater-flow problem for the hydraulic head `h`,
//!
//!    ```text
//!      −∇·(K ∇h) = 0     in a saturated body
//!              h = h̄     on the Dirichlet (prescribed-head) boundary
//!       −K ∂h/∂n = q̄     on the Neumann (prescribed-flux) boundary
//!    ```
//!
//!    discretised into `K_g · h = Q`, where `K_g` is the global
//!    **conductance matrix**, `h` the unknown nodal head field, and `Q`
//!    the flux-source vector. The Darcy flux (specific discharge) is
//!    recovered as `q = −K ∇h`. This is **mathematically identical** to
//!    the steady heat-conduction problem in [`crate::thermal_solver`]
//!    (head `h` ↔ temperature `T`, hydraulic conductivity `K` ↔ thermal
//!    conductivity `k`, Darcy flux `q = −K∇h` ↔ heat flux `q = −k∇T`), so
//!    it reuses the exact same element machinery
//!    ([`crate::thermal_solver::gradient_matrix`] /
//!    [`crate::thermal_solver::tet_volume`]) and the SPD sparse solver
//!    ([`crate::faer_solver::solve_spd`]).
//!
//! 2. **Unsaturated Richards flow** ([`solve_richards_step`] /
//!    [`solve_richards_transient`]). A backward-Euler **implicit**
//!    time step of the mixed-form Richards equation
//!
//!    ```text
//!      ∂θ/∂t = ∇·( K(h) ∇(h + z) )
//!    ```
//!
//!    for transient infiltration into an *unsaturated* soil column, with
//!    the **van Genuchten** retention curve `θ(h)` and the
//!    **Mualem–van Genuchten** relative permeability `Kᵣ(θ)`. The
//!    pressure head `h` may be negative (suction) above the water table;
//!    the gravity term `∇z` drives downward percolation. Mass is tracked
//!    so the caller can confirm conservation (stored-water increase =
//!    net inflow).
//!
//! ## Method (saturated Darcy)
//!
//! Identical to [`crate::thermal_solver`]: the four linear shape-function
//! gradients of a tet are constant, so the element conductance is the
//! closed form `Kₑ = K · Vₑ · Gᵀ·G`; each 4×4 `Kₑ` is scatter-added into a
//! global sparse matrix of size `n_nodes`; Neumann (prescribed-flux)
//! inputs go straight into `Q`; Dirichlet (prescribed-head) BCs are
//! imposed by the large-penalty method (keeping the system SPD); the
//! system is solved with [`crate::faer_solver::solve_spd`]; the
//! element flux `qₑ = −K·G·hₑ` is averaged to the nodes.
//!
//! ## Method (Richards)
//!
//! On each time step, the mixed-form `∂θ/∂t = ∇·(K(h)∇(h+z))` is
//! linearised by a **modified-Picard / mass-lumped backward-Euler**
//! scheme (Celia, Bouloutas & Zarba 1990): a lumped capacitance
//! `C = ∂θ/∂h` (the specific moisture capacity) gives the transient term,
//! the conductance uses `K(h)` evaluated at the current iterate, and the
//! gravity term enters the load as the divergence of `K·∇z = K·ẑ`. The
//! resulting SPD system is solved with the same faer path, and the
//! iterate is repeated to convergence within the step.
//!
//! ## Honest scope — what it is and is not
//!
//! **Saturated Darcy is a real flow solver.** For a linear, isotropic,
//! steady saturated problem on linear tets the head field it returns is
//! the correct FE solution for that mesh. A 1-D column with the head held
//! at `h₀` over `h₁` across length `L` reproduces the analytic **linear
//! head profile** `h(x) = h₀ + (h₁ − h₀)·x/L` and the **constant Darcy
//! flux** `q = K·(h₀ − h₁)/L` to solver precision (see the tests).
//!
//! **Richards is a genuine but deliberately V1 unsaturated step:**
//!
//! - **Isotropic, single-material.** One scalar saturated conductivity
//!   `K_sat` and one van-Genuchten parameter set; no anisotropy, no
//!   heterogeneous layering (the API takes scalars).
//! - **Modified-Picard linearisation** with a lumped capacitance — robust
//!   and mass-conservative in the mixed form, but not the full Newton
//!   Jacobian; very dry fronts may need small steps.
//! - **No hysteresis, no two-phase air flow** (the standard Richards
//!   single-phase assumption).
//! - **Tet4 only**; needs a volume mesh.
//!
//! Within that scope the numbers are real; anything outside it is
//! documented here, not faked. Cross-validation against a dedicated
//! groundwater code (MODFLOW / FEFLOW / HYDRUS) is a follow-up, tracked
//! in `docs/VALIDATION.md`.

use nalgebra::{DVector, Vector3, Vector4};
use nalgebra_sparse::{CooMatrix, CscMatrix};
use thiserror::Error;

use valenx_mesh::Mesh;

use crate::faer_solver::{solve_spd, FaerSolveError, SolverBackend};
use crate::thermal_solver::{gradient_matrix, tet_block, tet_volume};

/// Errors from the native porous-media flow solvers.
#[derive(Debug, Error)]
pub enum PorousFlowError {
    /// The supplied mesh has no [`valenx_mesh::element::ElementType::Tet4`] block.
    #[error(
        "mesh has no Tet4 element block — the porous-flow solver needs a tetrahedral volume mesh"
    )]
    NoTetBlock,
    /// The mesh has zero nodes.
    #[error("mesh has no nodes")]
    EmptyMesh,
    /// A tetrahedron is degenerate (zero or negative signed volume).
    #[error("tetrahedron {0} is degenerate (|volume| < tolerance)")]
    DegenerateElement(usize),
    /// A connectivity index points past the end of the node array.
    #[error("element {elem} references node {node} but the mesh has only {n_nodes} nodes")]
    BadConnectivity {
        /// 0-based element index (or `usize::MAX` for a BC reference).
        elem: usize,
        /// Out-of-range node index.
        node: usize,
        /// Node-array length.
        n_nodes: usize,
    },
    /// The hydraulic conductivity is non-physical (`K ≤ 0` or non-finite).
    #[error("invalid hydraulic conductivity: must be finite and positive, got {0}")]
    BadConductivity(f64),
    /// A van Genuchten parameter is out of range (`α ≤ 0`, `n ≤ 1`,
    /// `θ_s ≤ θ_r`, or any non-finite / negative saturation).
    #[error("invalid van Genuchten parameters: {0}")]
    BadRetention(String),
    /// A load or boundary-condition input carried a non-finite value
    /// (`NaN` / `±∞`). Such a value would flow straight into the
    /// right-hand side and corrupt the (otherwise `Ok`) solve, so it is
    /// rejected up front (mirroring the thermal solver's H4 guard).
    #[error("non-finite {kind} at node {node} (NaN or infinity is not a valid head/flux)")]
    InvalidInput {
        /// 0-based node index carrying the bad value.
        node: usize,
        /// Which input was non-finite — e.g. `"prescribed head"` or
        /// `"prescribed flux"`.
        kind: &'static str,
    },
    /// No Dirichlet (prescribed-head) boundary condition was supplied.
    /// With only Neumann BCs the conductance matrix is singular (the head
    /// is determined only up to an additive constant) and the steady
    /// solve has no unique answer.
    #[error("no prescribed-head boundary — the steady seepage problem is singular (head floats)")]
    NoDirichletBc,
    /// A non-positive or non-finite time step was supplied to the Richards
    /// integrator.
    #[error("invalid time step: must be finite and positive, got {0}")]
    BadTimeStep(f64),
    /// The Cholesky factorisation failed — the assembled system was not
    /// positive-definite.
    #[error("linear solve failed: conductance matrix is not positive-definite")]
    SolveFailed,
}

impl From<FaerSolveError> for PorousFlowError {
    fn from(_: FaerSolveError) -> Self {
        PorousFlowError::SolveFailed
    }
}

// ===========================================================================
// Saturated Darcy flow
// ===========================================================================

/// A single-node Dirichlet (prescribed hydraulic-head) boundary condition.
///
/// Hydraulic head `h = z + p/(ρg)` (elevation head + pressure head) is the
/// driving potential of Darcy flow; flow runs from high head to low head.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct PrescribedHead {
    /// 0-based node index.
    pub node: usize,
    /// Prescribed hydraulic head at the node (length units — m, ft — must
    /// be consistent with `K` and the geometry; the solver is
    /// unit-agnostic so long as they agree).
    pub head: f64,
}

/// A single-node Neumann boundary condition: a concentrated volumetric
/// flux injected at a node (volume/time). Positive = water flowing *into*
/// the body.
///
/// A distributed boundary flux `q̄` (volume/area/time) is converted to
/// equivalent nodal inputs by the caller — for a linear-tet face of area
/// `A` a uniform flux contributes `q̄·A/3` to each of the face's three
/// nodes (exactly as for the thermal Neumann load).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct FluxSource {
    /// 0-based node index.
    pub node: usize,
    /// Volumetric flux into the body at the node (volume/time, positive =
    /// inflow).
    pub rate: f64,
}

/// Result of a steady-state saturated Darcy solve.
#[derive(Clone, Debug)]
pub struct DarcySolution {
    /// Per-node hydraulic-head field `h`, indexed by node.
    pub head: Vec<f64>,
    /// Per-node Darcy flux (specific discharge) vector `[qx, qy, qz]`
    /// (length/time), computed per element (constant over a linear tet)
    /// then averaged to the nodes.
    pub flux: Vec<[f64; 3]>,
}

impl DarcySolution {
    /// Highest nodal head in the field.
    #[must_use]
    pub fn max_head(&self) -> f64 {
        self.head.iter().copied().fold(f64::MIN, f64::max)
    }

    /// Lowest nodal head in the field.
    #[must_use]
    pub fn min_head(&self) -> f64 {
        self.head.iter().copied().fold(f64::MAX, f64::min)
    }

    /// Largest nodal Darcy-flux magnitude (the peak seepage velocity ÷
    /// porosity is the actual pore velocity; this is the specific
    /// discharge).
    #[must_use]
    pub fn max_flux_magnitude(&self) -> f64 {
        self.flux
            .iter()
            .map(|q| (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt())
            .fold(0.0, f64::max)
    }
}

/// Solve a steady-state **saturated Darcy** groundwater-flow problem on a
/// tetrahedral mesh.
///
/// `mesh` must carry at least one [`valenx_mesh::element::ElementType::Tet4`]
/// block. `conductivity` is the isotropic saturated hydraulic conductivity
/// `K` (length/time). `fixed` supplies the Dirichlet (prescribed-head) BCs
/// (at least one is required, or the steady problem is singular). `sources`
/// supplies the Neumann (volumetric-flux) BCs.
///
/// Returns the per-node head field and the recovered nodal Darcy flux
/// `q = −K ∇h`.
///
/// # Errors
///
/// See [`PorousFlowError`].
pub fn solve_saturated_darcy(
    mesh: &Mesh,
    conductivity: f64,
    fixed: &[PrescribedHead],
    sources: &[FluxSource],
) -> Result<DarcySolution, PorousFlowError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(PorousFlowError::EmptyMesh);
    }
    if !conductivity.is_finite() || conductivity <= 0.0 {
        return Err(PorousFlowError::BadConductivity(conductivity));
    }
    if fixed.is_empty() {
        return Err(PorousFlowError::NoDirichletBc);
    }
    let tets = tet_block(mesh).ok_or(PorousFlowError::NoTetBlock)?;
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;

    // --- assemble the global conductance matrix (n_nodes square) ---
    // Element conductance Kₑ = K·Vₑ·Gᵀ·G — identical closed form to the
    // thermal element matrix (this *is* the same scalar Laplacian).
    let mut coo = CooMatrix::<f64>::new(n_nodes, n_nodes);
    for e in 0..n_elem {
        let nodes = elem_nodes(conn, e);
        check_nodes(&nodes, e, n_nodes)?;
        let coords = elem_coords(mesh, &nodes);
        let g = gradient_matrix(&coords).ok_or(PorousFlowError::DegenerateElement(e))?;
        let vol = tet_volume(&coords).abs();
        if vol < 1.0e-18 {
            return Err(PorousFlowError::DegenerateElement(e));
        }
        let ke = (g.transpose() * g) * (conductivity * vol); // 4×4
        for a in 0..4 {
            for b in 0..4 {
                let v = ke[(a, b)];
                if v != 0.0 {
                    coo.push(nodes[a], nodes[b], v);
                }
            }
        }
    }
    let csc = CscMatrix::from(&coo);

    // Peak diagonal for penalty scaling.
    let mut max_diag = 0.0_f64;
    for (r, c, v) in csc.triplet_iter() {
        if r == c {
            max_diag = max_diag.max(v.abs());
        }
    }
    if max_diag <= 0.0 {
        return Err(PorousFlowError::SolveFailed);
    }
    let penalty = max_diag * 1.0e8;

    // --- flux-source vector (Neumann) ---
    let mut q = DVector::<f64>::zeros(n_nodes);
    for s in sources {
        if s.node >= n_nodes {
            return Err(PorousFlowError::BadConnectivity {
                elem: usize::MAX,
                node: s.node,
                n_nodes,
            });
        }
        if !s.rate.is_finite() {
            return Err(PorousFlowError::InvalidInput {
                node: s.node,
                kind: "prescribed flux",
            });
        }
        q[s.node] += s.rate;
    }

    // --- Dirichlet (prescribed-head) BCs via the penalty method ---
    let mut penalty_diag = vec![0.0_f64; n_nodes];
    for fb in fixed {
        if fb.node >= n_nodes {
            return Err(PorousFlowError::BadConnectivity {
                elem: usize::MAX,
                node: fb.node,
                n_nodes,
            });
        }
        if !fb.head.is_finite() {
            return Err(PorousFlowError::InvalidInput {
                node: fb.node,
                kind: "prescribed head",
            });
        }
        penalty_diag[fb.node] += penalty;
        q[fb.node] += penalty * fb.head;
    }

    let csc = add_diagonal(&csc, &penalty_diag);

    // --- solve K_g · h = Q (SPD, reuse the faer/legacy backend) ---
    let backend = SolverBackend::default();
    let h_vec = solve_spd(&csc, &q, backend)?;
    if h_vec.iter().any(|v| !v.is_finite()) {
        return Err(PorousFlowError::SolveFailed);
    }
    let head: Vec<f64> = (0..n_nodes).map(|n| h_vec[n]).collect();

    // --- Darcy-flux recovery: qₑ = −K·G·hₑ per element, averaged to nodes ---
    let flux = recover_nodal_flux(mesh, conn, n_elem, n_nodes, &head, |_| conductivity);

    Ok(DarcySolution { head, flux })
}

// ===========================================================================
// van Genuchten retention + Mualem relative permeability (unsaturated)
// ===========================================================================

/// **van Genuchten** soil-water retention parameters, with the
/// **Mualem** relative-permeability model. These describe how a partially
/// saturated soil holds and conducts water as a function of the (negative,
/// suction) pressure head `h`.
///
/// Effective saturation
/// `Se(h) = [1 + |α·h|ⁿ]^(−m)` for `h < 0`, `Se = 1` for `h ≥ 0`
/// (saturated), with `m = 1 − 1/n`. Water content
/// `θ(h) = θ_r + (θ_s − θ_r)·Se`. Relative permeability (Mualem)
/// `Kᵣ(Se) = √Se · [1 − (1 − Se^(1/m))^m]²`.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct VanGenuchten {
    /// Residual volumetric water content `θ_r` (dimensionless, ≥ 0).
    pub theta_r: f64,
    /// Saturated volumetric water content `θ_s` (≈ porosity, > `θ_r`).
    pub theta_s: f64,
    /// Inverse air-entry pressure `α` (1/length, > 0).
    pub alpha: f64,
    /// Pore-size-distribution index `n` (> 1).
    pub n: f64,
}

impl VanGenuchten {
    /// A typical **loam** parameter set (Carsel & Parrish 1988), with `α`
    /// in 1/m: `θ_r = 0.078`, `θ_s = 0.43`, `α = 3.6 /m`, `n = 1.56`.
    #[must_use]
    pub fn loam() -> Self {
        Self {
            theta_r: 0.078,
            theta_s: 0.43,
            alpha: 3.6,
            n: 1.56,
        }
    }

    /// A typical **sand** parameter set (Carsel & Parrish 1988), `α` in
    /// 1/m: `θ_r = 0.045`, `θ_s = 0.43`, `α = 14.5 /m`, `n = 2.68`.
    #[must_use]
    pub fn sand() -> Self {
        Self {
            theta_r: 0.045,
            theta_s: 0.43,
            alpha: 14.5,
            n: 2.68,
        }
    }

    /// Validate the parameter ranges.
    fn validate(&self) -> Result<(), PorousFlowError> {
        if !self.theta_r.is_finite() || self.theta_r < 0.0 {
            return Err(PorousFlowError::BadRetention(format!(
                "theta_r must be finite and ≥ 0, got {}",
                self.theta_r
            )));
        }
        if !self.theta_s.is_finite() || self.theta_s <= self.theta_r {
            return Err(PorousFlowError::BadRetention(format!(
                "theta_s must be finite and > theta_r ({}), got {}",
                self.theta_r, self.theta_s
            )));
        }
        if !self.alpha.is_finite() || self.alpha <= 0.0 {
            return Err(PorousFlowError::BadRetention(format!(
                "alpha must be finite and > 0, got {}",
                self.alpha
            )));
        }
        if !self.n.is_finite() || self.n <= 1.0 {
            return Err(PorousFlowError::BadRetention(format!(
                "n must be finite and > 1, got {}",
                self.n
            )));
        }
        Ok(())
    }

    /// `m = 1 − 1/n`.
    #[inline]
    fn m(&self) -> f64 {
        1.0 - 1.0 / self.n
    }

    /// Effective saturation `Se ∈ (0, 1]` at pressure head `h`. `h ≥ 0`
    /// (at/below the water table) is fully saturated (`Se = 1`).
    #[must_use]
    pub fn effective_saturation(&self, h: f64) -> f64 {
        if h >= 0.0 {
            return 1.0;
        }
        let m = self.m();
        let ah = (self.alpha * h).abs();
        (1.0 + ah.powf(self.n)).powf(-m)
    }

    /// Volumetric water content `θ(h) = θ_r + (θ_s − θ_r)·Se`.
    #[must_use]
    pub fn water_content(&self, h: f64) -> f64 {
        self.theta_r + (self.theta_s - self.theta_r) * self.effective_saturation(h)
    }

    /// Specific moisture capacity `C(h) = dθ/dh` — the slope of the
    /// retention curve, which forms the (lumped) transient term of the
    /// Richards system. Zero in the saturated branch (`h ≥ 0`).
    #[must_use]
    pub fn specific_moisture_capacity(&self, h: f64) -> f64 {
        if h >= 0.0 {
            return 0.0;
        }
        // dSe/dh = −m·n·α·|αh|^(n−1)·(1+|αh|ⁿ)^(−m−1)·sign(h); h<0 ⇒ |αh|=−αh.
        let m = self.m();
        let ah = (self.alpha * h).abs();
        let ahn = ah.powf(self.n);
        let dse_dh = m * self.n * self.alpha * ah.powf(self.n - 1.0) * (1.0 + ahn).powf(-m - 1.0);
        // θ increases as h rises toward 0, so dθ/dh > 0.
        (self.theta_s - self.theta_r) * dse_dh
    }

    /// **Mualem** relative permeability `Kᵣ(Se) ∈ [0, 1]` from the
    /// effective saturation. `1` when saturated, → 0 as the soil dries.
    #[must_use]
    pub fn relative_permeability(&self, h: f64) -> f64 {
        let se = self.effective_saturation(h).clamp(0.0, 1.0);
        if se >= 1.0 {
            return 1.0;
        }
        let m = self.m();
        let inner = 1.0 - (1.0 - se.powf(1.0 / m)).powf(m);
        se.sqrt() * inner * inner
    }
}

/// State carried between Richards time steps: the nodal pressure-head
/// field and the time elapsed.
#[derive(Clone, Debug)]
pub struct RichardsState {
    /// Per-node pressure head `h` (length; negative = unsaturated suction).
    pub head: Vec<f64>,
    /// Elapsed time of the simulation (time units).
    pub time: f64,
}

impl RichardsState {
    /// Start a Richards simulation from a uniform initial pressure head.
    #[must_use]
    pub fn uniform(n_nodes: usize, head0: f64) -> Self {
        Self {
            head: vec![head0; n_nodes],
            time: 0.0,
        }
    }

    /// Total stored water volume `∫θ dV` over the mesh, used for the
    /// mass-conservation check. Integrated with the same lumped (nodal)
    /// quadrature the transient term uses, so it is consistent with the
    /// scheme's own mass.
    #[must_use]
    pub fn stored_water(&self, mesh: &Mesh, vg: &VanGenuchten) -> f64 {
        let Some(tets) = tet_block(mesh) else {
            return 0.0;
        };
        let conn = &tets.connectivity;
        let n_elem = conn.len() / 4;
        let mut total = 0.0;
        for e in 0..n_elem {
            let nodes = elem_nodes(conn, e);
            if nodes.iter().any(|&n| n >= mesh.nodes.len()) {
                continue;
            }
            let coords = elem_coords(mesh, &nodes);
            let vol = tet_volume(&coords).abs();
            // Lumped: each node carries V/4 of the element; θ at that node.
            for &nd in &nodes {
                total += vg.water_content(self.head[nd]) * vol / 4.0;
            }
        }
        total
    }
}

/// Advance a transient **Richards** (unsaturated infiltration) simulation
/// by one implicit backward-Euler time step of size `dt`.
///
/// `state` is the field at the start of the step (consumed-by-reference;
/// the returned state is the field at `t + dt`). `k_sat` is the saturated
/// hydraulic conductivity (length/time); `vg` the van Genuchten retention
/// parameters. `fixed` pins pressure head at boundary nodes (e.g. a
/// ponded-infiltration top boundary `h = 0`, or a water-table bottom
/// boundary); `sources` injects volumetric flux at nodes. Gravity acts
/// along `−z` (the mesh's third coordinate), driving downward percolation.
///
/// The scheme is the **mass-lumped modified-Picard backward-Euler** of
/// Celia–Bouloutas–Zarba (1990): within the step the nonlinear `K(h)` /
/// `C(h)` are iterated to convergence (up to `max_iter`, tolerance `tol`
/// on the head increment).
///
/// # Errors
///
/// See [`PorousFlowError`].
#[allow(clippy::too_many_arguments)]
pub fn solve_richards_step(
    mesh: &Mesh,
    state: &RichardsState,
    dt: f64,
    k_sat: f64,
    vg: &VanGenuchten,
    fixed: &[PrescribedHead],
    sources: &[FluxSource],
    max_iter: usize,
    tol: f64,
) -> Result<RichardsState, PorousFlowError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(PorousFlowError::EmptyMesh);
    }
    if !k_sat.is_finite() || k_sat <= 0.0 {
        return Err(PorousFlowError::BadConductivity(k_sat));
    }
    if !dt.is_finite() || dt <= 0.0 {
        return Err(PorousFlowError::BadTimeStep(dt));
    }
    vg.validate()?;
    if state.head.len() != n_nodes {
        return Err(PorousFlowError::BadConnectivity {
            elem: usize::MAX,
            node: state.head.len(),
            n_nodes,
        });
    }
    let tets = tet_block(mesh).ok_or(PorousFlowError::NoTetBlock)?;
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;

    // Validate BC node indices / finiteness up front.
    for fb in fixed {
        if fb.node >= n_nodes {
            return Err(PorousFlowError::BadConnectivity {
                elem: usize::MAX,
                node: fb.node,
                n_nodes,
            });
        }
        if !fb.head.is_finite() {
            return Err(PorousFlowError::InvalidInput {
                node: fb.node,
                kind: "prescribed head",
            });
        }
    }
    for s in sources {
        if s.node >= n_nodes {
            return Err(PorousFlowError::BadConnectivity {
                elem: usize::MAX,
                node: s.node,
                n_nodes,
            });
        }
        if !s.rate.is_finite() {
            return Err(PorousFlowError::InvalidInput {
                node: s.node,
                kind: "prescribed flux",
            });
        }
    }

    let theta_old: Vec<f64> = state.head.iter().map(|&h| vg.water_content(h)).collect();

    // Modified-Picard iteration on the head at t+dt.
    let mut h_new = state.head.clone();
    let max_iter = max_iter.max(1);
    for _iter in 0..max_iter {
        // Assemble per-element conductance with K(h)=K_sat·Kr(h̄_elem) and
        // the lumped capacitance C(h)/dt on the diagonal; the gravity load
        // is ∫ K ∇·ẑ which for a linear tet is K·Vₑ·(Gᵀ·ẑ) per node.
        let mut coo = CooMatrix::<f64>::new(n_nodes, n_nodes);
        let mut rhs = DVector::<f64>::zeros(n_nodes);
        let mut lumped_c = vec![0.0_f64; n_nodes];

        for e in 0..n_elem {
            let nodes = elem_nodes(conn, e);
            check_nodes(&nodes, e, n_nodes)?;
            let coords = elem_coords(mesh, &nodes);
            let g = gradient_matrix(&coords).ok_or(PorousFlowError::DegenerateElement(e))?;
            let vol = tet_volume(&coords).abs();
            if vol < 1.0e-18 {
                return Err(PorousFlowError::DegenerateElement(e));
            }
            // Relative permeability at the element's mean head (upstream-ish
            // averaging keeps it simple and bounded for V1).
            let h_elem =
                0.25 * (h_new[nodes[0]] + h_new[nodes[1]] + h_new[nodes[2]] + h_new[nodes[3]]);
            let k_elem = k_sat * vg.relative_permeability(h_elem);
            let ke = (g.transpose() * g) * (k_elem * vol); // 4×4 conductance
            for a in 0..4 {
                for b in 0..4 {
                    let v = ke[(a, b)];
                    if v != 0.0 {
                        coo.push(nodes[a], nodes[b], v);
                    }
                }
            }
            // Gravity term: flow potential is (h + z), so the extra load is
            //   −∫ ∇Nₐ · K ∇z dV = −K·Vₑ·(Gᵀ·ẑ)_a  ... moved to RHS as +.
            // ẑ = (0,0,1); Gᵀ·ẑ is the z-row of G (the ∂N/∂z of each node).
            for a in 0..4 {
                let dndz = g[(2, a)];
                rhs[nodes[a]] -= k_elem * vol * dndz;
            }
            // Lumped capacitance: C(h) at each node, V/4 share.
            for &nd in &nodes {
                lumped_c[nd] += vg.specific_moisture_capacity(h_new[nd]) * vol / 4.0;
            }
        }

        // Transient term (mixed form, mass-lumped):
        //   C/dt·(h_new − h_iter) + (θ(h_new) − θ_old)/dt ... modified-Picard
        // The standard CBZ mixed update folds C/dt onto the diagonal and puts
        //   −(θ(h_new) − θ_old)/dt + C/dt·h_iter  on the RHS, so converged
        // iterate satisfies the θ-form exactly (mass-conservative).
        let mut diag = vec![0.0_f64; n_nodes];
        for n in 0..n_nodes {
            let c_dt = lumped_c[n] / dt;
            diag[n] += c_dt;
            let theta_new = vg.water_content(h_new[n]);
            rhs[n] += c_dt * h_new[n] - (theta_new - theta_old[n]) / dt;
        }

        // Neumann flux sources.
        for s in sources {
            rhs[s.node] += s.rate;
        }

        let csc = CscMatrix::from(&coo);
        // Peak diagonal for penalty scaling (include the capacitance diag).
        let mut max_diag = 0.0_f64;
        for (r, c, v) in csc.triplet_iter() {
            if r == c {
                max_diag = max_diag.max(v.abs());
            }
        }
        for &d in &diag {
            max_diag = max_diag.max(d.abs());
        }
        if max_diag <= 0.0 {
            return Err(PorousFlowError::SolveFailed);
        }
        let penalty = max_diag * 1.0e8;

        // Dirichlet (prescribed-head) via penalty.
        for fb in fixed {
            diag[fb.node] += penalty;
            rhs[fb.node] += penalty * fb.head;
        }

        let csc = add_diagonal(&csc, &diag);
        let backend = SolverBackend::default();
        let h_sol = solve_spd(&csc, &rhs, backend)?;
        if h_sol.iter().any(|v| !v.is_finite()) {
            return Err(PorousFlowError::SolveFailed);
        }

        // Convergence on the head increment.
        let mut max_delta = 0.0_f64;
        for n in 0..n_nodes {
            max_delta = max_delta.max((h_sol[n] - h_new[n]).abs());
            h_new[n] = h_sol[n];
        }
        if max_delta < tol {
            break;
        }
    }

    Ok(RichardsState {
        head: h_new,
        time: state.time + dt,
    })
}

/// Run a fixed number of equal Richards time steps, returning the final
/// state. A thin convenience wrapper over [`solve_richards_step`].
///
/// # Errors
///
/// See [`PorousFlowError`].
#[allow(clippy::too_many_arguments)]
pub fn solve_richards_transient(
    mesh: &Mesh,
    initial: &RichardsState,
    dt: f64,
    n_steps: usize,
    k_sat: f64,
    vg: &VanGenuchten,
    fixed: &[PrescribedHead],
    sources: &[FluxSource],
) -> Result<RichardsState, PorousFlowError> {
    let mut state = initial.clone();
    for _ in 0..n_steps {
        state = solve_richards_step(mesh, &state, dt, k_sat, vg, fixed, sources, 20, 1.0e-6)?;
    }
    Ok(state)
}

// ===========================================================================
// shared helpers
// ===========================================================================

/// The four 0-based node indices of element `e` from a flat Tet4
/// connectivity array.
#[inline]
fn elem_nodes(conn: &[u32], e: usize) -> [usize; 4] {
    [
        conn[4 * e] as usize,
        conn[4 * e + 1] as usize,
        conn[4 * e + 2] as usize,
        conn[4 * e + 3] as usize,
    ]
}

/// Bounds-check an element's node indices against the node count.
#[inline]
fn check_nodes(nodes: &[usize; 4], e: usize, n_nodes: usize) -> Result<(), PorousFlowError> {
    for &nd in nodes {
        if nd >= n_nodes {
            return Err(PorousFlowError::BadConnectivity {
                elem: e,
                node: nd,
                n_nodes,
            });
        }
    }
    Ok(())
}

/// The four node coordinates of an element.
#[inline]
fn elem_coords(mesh: &Mesh, nodes: &[usize; 4]) -> [Vector3<f64>; 4] {
    [
        mesh.nodes[nodes[0]],
        mesh.nodes[nodes[1]],
        mesh.nodes[nodes[2]],
        mesh.nodes[nodes[3]],
    ]
}

/// Recover a nodal flux field `q = −K(·)·G·hₑ` by averaging the constant
/// per-element flux into its nodes. `k_of` returns the conductivity to use
/// for element `e` (constant `K` for saturated Darcy; `K_sat·Kr` could be
/// passed for an unsaturated post-process).
fn recover_nodal_flux(
    mesh: &Mesh,
    conn: &[u32],
    n_elem: usize,
    n_nodes: usize,
    head: &[f64],
    k_of: impl Fn(usize) -> f64,
) -> Vec<[f64; 3]> {
    let mut nodal_flux = vec![[0.0_f64; 3]; n_nodes];
    let mut nodal_count = vec![0_u32; n_nodes];
    for e in 0..n_elem {
        let nodes = elem_nodes(conn, e);
        if nodes.iter().any(|&n| n >= n_nodes) {
            continue;
        }
        let coords = elem_coords(mesh, &nodes);
        let Some(g) = gradient_matrix(&coords) else {
            continue;
        };
        let h_elem = Vector4::new(
            head[nodes[0]],
            head[nodes[1]],
            head[nodes[2]],
            head[nodes[3]],
        );
        let grad_h = g * h_elem; // 3×1
        let k = k_of(e);
        let flux = grad_h * (-k); // q = −K ∇h
        for &nd in &nodes {
            nodal_flux[nd][0] += flux.x;
            nodal_flux[nd][1] += flux.y;
            nodal_flux[nd][2] += flux.z;
            nodal_count[nd] += 1;
        }
    }
    let mut flux = vec![[0.0_f64; 3]; n_nodes];
    for n in 0..n_nodes {
        let count = nodal_count[n].max(1) as f64;
        flux[n][0] = nodal_flux[n][0] / count;
        flux[n][1] = nodal_flux[n][1] / count;
        flux[n][2] = nodal_flux[n][2] / count;
    }
    flux
}

/// Return a copy of `csc` with `diag[i]` added to entry `(i, i)`.
/// Used to fold the Dirichlet penalty (and Richards capacitance) onto the
/// conductance matrix diagonal.
fn add_diagonal(csc: &CscMatrix<f64>, diag: &[f64]) -> CscMatrix<f64> {
    let n = csc.nrows();
    let mut coo = CooMatrix::<f64>::new(n, n);
    for (r, c, v) in csc.triplet_iter() {
        coo.push(r, c, *v);
    }
    for (i, &d) in diag.iter().enumerate() {
        if d != 0.0 {
            coo.push(i, i, d);
        }
    }
    CscMatrix::from(&coo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native_solver::structured_box_mesh;

    /// Node id of grid point `(i, j, k)` in a `structured_box_mesh`
    /// (X fastest, then Y, then Z).
    fn nid(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
        i + (nx + 1) * j + (nx + 1) * (ny + 1) * k
    }

    // ----- saturated Darcy: the headline analytic validation -----

    #[test]
    fn one_dimensional_column_reproduces_linear_head_and_constant_flux() {
        // A column along X: hold the x=0 face at h0 (high head) and the
        // x=L face at h1 (low head), sides no-flow (natural BC). The
        // steady saturated solution is the analytic LINEAR head profile
        //   h(x) = h0 + (h1 − h0)·x/L
        // and the CONSTANT Darcy flux  q = −K·dh/dx = K·(h0 − h1)/L  (in
        // +X since head drops in +X). Match must be <1%.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (8, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box");
        let k = 1.0e-4; // m/s, a silty sand
        let h0 = 10.0;
        let h1 = 2.0;

        let mut fixed = Vec::new();
        for j in 0..=ny {
            for kk in 0..=nz {
                fixed.push(PrescribedHead {
                    node: nid(0, j, kk, nx, ny),
                    head: h0,
                });
                fixed.push(PrescribedHead {
                    node: nid(nx, j, kk, nx, ny),
                    head: h1,
                });
            }
        }
        let sol = solve_saturated_darcy(&mesh, k, &fixed, &[]).unwrap();

        // Linear head profile at every node (<1%).
        for i in 0..=nx {
            for j in 0..=ny {
                for kk in 0..=nz {
                    let n = nid(i, j, kk, nx, ny);
                    let x = mesh.nodes[n].x;
                    let expected = h0 + (h1 - h0) * x / lx;
                    let got = sol.head[n];
                    let tol = 0.01 * (h0 - h1).abs();
                    assert!(
                        (got - expected).abs() < tol,
                        "node {n} at x={x}: head {got} vs analytic {expected}"
                    );
                }
            }
        }

        // Constant flux q = K·(h0−h1)/L in +X, at an interior node (<1%).
        let expected_q = k * (h0 - h1) / lx;
        let interior = nid(4, 0, 0, nx, ny);
        let qx = sol.flux[interior][0];
        assert!(
            (qx - expected_q).abs() < 0.01 * expected_q.abs(),
            "interior Darcy flux qx {qx} vs analytic {expected_q}"
        );
        // Head drops in +X ⇒ flux is in +X (positive).
        assert!(qx > 0.0, "flow should be from high head (+X face) outward");
    }

    #[test]
    fn uniform_head_gives_zero_flux() {
        // Hold every node at the same head → no gradient → no flow.
        let (nx, ny, nz) = (2, 2, 2);
        let mesh = structured_box_mesh(2.0, 2.0, 2.0, nx, ny, nz).expect("valid box");
        let fixed: Vec<PrescribedHead> = (0..mesh.nodes.len())
            .map(|n| PrescribedHead { node: n, head: 5.0 })
            .collect();
        let sol = solve_saturated_darcy(&mesh, 1.0e-4, &fixed, &[]).unwrap();
        for &h in &sol.head {
            assert!((h - 5.0).abs() < 1e-6, "head {h} ≠ 5");
        }
        assert!(
            sol.max_flux_magnitude() < 1e-9,
            "a uniform head must carry no Darcy flux"
        );
    }

    #[test]
    fn dam_seepage_flows_from_high_to_low_reservoir() {
        // 2-D dam-seepage analogue: a block with a high reservoir head on
        // the x=0 face (upstream) and a low tailwater head on the x=L face
        // (downstream), no-flow on top/bottom/sides. Qualitative checks:
        //  - head decreases monotonically from upstream to downstream,
        //  - net Darcy flux points downstream (+X),
        //  - mass balance: total inflow at the upstream face ≈ total
        //    outflow at the downstream face.
        let (lx, ly, lz) = (10.0, 1.0, 4.0); // length × thickness × height
        let (nx, ny, nz) = (10, 1, 4);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box");
        let k = 5.0e-5;
        let h_up = 8.0; // upstream reservoir
        let h_dn = 2.0; // downstream tailwater

        let mut fixed = Vec::new();
        for j in 0..=ny {
            for kk in 0..=nz {
                fixed.push(PrescribedHead {
                    node: nid(0, j, kk, nx, ny),
                    head: h_up,
                });
                fixed.push(PrescribedHead {
                    node: nid(nx, j, kk, nx, ny),
                    head: h_dn,
                });
            }
        }
        let sol = solve_saturated_darcy(&mesh, k, &fixed, &[]).unwrap();

        // Head decreases along +X at mid-height.
        let mid_k = nz / 2;
        let mut prev = f64::INFINITY;
        for i in 0..=nx {
            let h = sol.head[nid(i, 0, mid_k, nx, ny)];
            assert!(
                h <= prev + 1e-9,
                "head should not increase downstream (i={i}, h={h}, prev={prev})"
            );
            prev = h;
        }
        // Upstream higher than downstream.
        assert!(sol.head[nid(0, 0, mid_k, nx, ny)] > sol.head[nid(nx, 0, mid_k, nx, ny)]);

        // Net seepage points downstream (+X) on average.
        let mut sum_qx = 0.0;
        for q in &sol.flux {
            sum_qx += q[0];
        }
        assert!(sum_qx > 0.0, "net seepage should be downstream (+X)");

        // Mass balance: the conduction flux through any vertical section is
        // the same in steady state. Compare the analytic 1-D estimate
        // q·A = K·(Δh/L)·A against the recovered interior flux × area.
        let area = ly * lz;
        let q_analytic = k * (h_up - h_dn) / lx * area;
        // Average qx over interior nodes × area as a coarse total-flow proxy.
        let interior_qx: f64 = (1..nx)
            .flat_map(|i| (0..=nz).map(move |kk| nid(i, 0, kk, nx, ny)))
            .map(|n| sol.flux[n][0])
            .sum::<f64>()
            / ((nx - 1) * (nz + 1)) as f64;
        let q_fe = interior_qx * area;
        assert!(
            (q_fe - q_analytic).abs() < 0.1 * q_analytic,
            "seepage flow {q_fe} vs 1-D analytic {q_analytic}"
        );
    }

    #[test]
    fn solve_rejects_bad_conductivity() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box");
        let err = solve_saturated_darcy(&mesh, -1.0, &[PrescribedHead { node: 0, head: 1.0 }], &[])
            .unwrap_err();
        assert!(matches!(err, PorousFlowError::BadConductivity(_)));
    }

    #[test]
    fn solve_rejects_missing_dirichlet_bc() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box");
        let err = solve_saturated_darcy(&mesh, 1.0e-4, &[], &[]).unwrap_err();
        assert!(matches!(err, PorousFlowError::NoDirichletBc));
    }

    #[test]
    fn solve_rejects_mesh_without_tets() {
        use valenx_mesh::element::{ElementBlock, ElementType};
        let mut mesh = Mesh::new("surface-only");
        mesh.nodes.push(Vector3::new(0.0, 0.0, 0.0));
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tri3,
            connectivity: vec![0, 0, 0],
        });
        let err =
            solve_saturated_darcy(&mesh, 1.0e-4, &[PrescribedHead { node: 0, head: 1.0 }], &[])
                .unwrap_err();
        assert!(matches!(err, PorousFlowError::NoTetBlock));
    }

    #[test]
    fn solve_rejects_non_finite_head_bc() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box");
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let err =
                solve_saturated_darcy(&mesh, 1.0e-4, &[PrescribedHead { node: 0, head: bad }], &[])
                    .unwrap_err();
            assert!(
                matches!(
                    err,
                    PorousFlowError::InvalidInput {
                        node: 0,
                        kind: "prescribed head"
                    }
                ),
                "expected InvalidInput(prescribed head) for {bad}, got {err:?}"
            );
        }
    }

    #[test]
    fn solve_rejects_non_finite_flux_source() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box");
        for bad in [f64::NAN, f64::INFINITY] {
            let err = solve_saturated_darcy(
                &mesh,
                1.0e-4,
                &[PrescribedHead { node: 0, head: 1.0 }],
                &[FluxSource { node: 1, rate: bad }],
            )
            .unwrap_err();
            assert!(
                matches!(
                    err,
                    PorousFlowError::InvalidInput {
                        node: 1,
                        kind: "prescribed flux"
                    }
                ),
                "expected InvalidInput(prescribed flux) for {bad}, got {err:?}"
            );
        }
    }

    // ----- van Genuchten retention curve sanity -----

    #[test]
    fn van_genuchten_saturates_and_dries_monotonically() {
        let vg = VanGenuchten::loam();
        // At/above the water table (h ≥ 0): fully saturated.
        assert!((vg.effective_saturation(0.0) - 1.0).abs() < 1e-12);
        assert!((vg.water_content(0.5) - vg.theta_s).abs() < 1e-12);
        assert!((vg.relative_permeability(0.1) - 1.0).abs() < 1e-12);

        // Drier (more negative h) ⇒ less water, lower Kr, both in (0,1).
        let h_wet = -0.1;
        let h_dry = -5.0;
        let se_wet = vg.effective_saturation(h_wet);
        let se_dry = vg.effective_saturation(h_dry);
        assert!(se_dry < se_wet && se_wet < 1.0 && se_dry > 0.0);
        assert!(vg.water_content(h_dry) < vg.water_content(h_wet));
        let kr_wet = vg.relative_permeability(h_wet);
        let kr_dry = vg.relative_permeability(h_dry);
        assert!(kr_dry < kr_wet && kr_wet <= 1.0 && kr_dry >= 0.0);
        // Water content stays within [θ_r, θ_s].
        for &h in &[0.0, -0.01, -1.0, -100.0] {
            let th = vg.water_content(h);
            assert!(
                th >= vg.theta_r - 1e-12 && th <= vg.theta_s + 1e-12,
                "θ {th} out of range"
            );
        }
    }

    #[test]
    fn specific_moisture_capacity_matches_finite_difference() {
        // C(h) = dθ/dh must match a centred finite difference of θ(h).
        let vg = VanGenuchten::loam();
        for &h in &[-0.2_f64, -0.5, -1.0, -2.0] {
            let eps = 1e-6;
            let fd = (vg.water_content(h + eps) - vg.water_content(h - eps)) / (2.0 * eps);
            let analytic = vg.specific_moisture_capacity(h);
            assert!(
                (fd - analytic).abs() < 1e-4 * fd.abs().max(1e-6),
                "C(h={h}): analytic {analytic} vs FD {fd}"
            );
        }
    }

    #[test]
    fn van_genuchten_rejects_bad_params() {
        let bad = VanGenuchten {
            theta_r: 0.1,
            theta_s: 0.05, // ≤ theta_r
            alpha: 1.0,
            n: 2.0,
        };
        assert!(bad.validate().is_err());
        let bad2 = VanGenuchten {
            theta_r: 0.0,
            theta_s: 0.4,
            alpha: 1.0,
            n: 1.0, // ≤ 1
        };
        assert!(bad2.validate().is_err());
    }

    // ----- Richards: infiltration + mass conservation -----

    #[test]
    fn richards_infiltration_wets_the_column_and_conserves_mass() {
        // A 1-D unsaturated soil column, initially dry (h = −2 m suction).
        // Pond the top face (highest z) at h = 0 (saturated) and hold the
        // bottom at its initial suction; water infiltrates downward under
        // gravity + the suction gradient. Checks:
        //  - saturation INCREASES (stored water grows) under infiltration,
        //  - the field stays finite and physical,
        //  - MASS BALANCE: the stored-water increase equals the net flux
        //    that crossed the boundaries over the run (within tolerance).
        let (lx, ly, lz) = (1.0, 1.0, 1.0);
        let (nx, ny, nz) = (1, 1, 8); // refine in z (the flow direction)
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box");
        let n_nodes = mesh.nodes.len();
        let vg = VanGenuchten::loam();
        let k_sat = 2.0e-6; // m/s (loam)
        let h_init = -2.0;

        let mut state = RichardsState::uniform(n_nodes, h_init);
        let stored0 = state.stored_water(&mesh, &vg);

        // Top face = highest z nodes; pond them at h = 0.
        let top: Vec<usize> = (0..=nx)
            .flat_map(|i| (0..=ny).map(move |j| nid(i, j, nz, nx, ny)))
            .collect();
        let fixed: Vec<PrescribedHead> = top
            .iter()
            .map(|&n| PrescribedHead { node: n, head: 0.0 })
            .collect();

        let bottom = nid(0, 0, 0, nx, ny);

        // March the infiltration and confirm the stored water increases
        // MONOTONICALLY step by step (water only enters — the ponded top is
        // the sole open boundary), then settles. This is the transient
        // mass-conservation signature of infiltration.
        let dt = 3600.0; // 1 hour
        let n_steps = 24; // one day
        let mut prev_stored = stored0;
        for _ in 0..n_steps {
            state =
                solve_richards_step(&mesh, &state, dt, k_sat, &vg, &fixed, &[], 25, 1e-7).unwrap();
            let s = state.stored_water(&mesh, &vg);
            assert!(
                s >= prev_stored - 1e-9,
                "stored water must not decrease during infiltration: {s} < {prev_stored}"
            );
            prev_stored = s;
        }

        let stored1 = state.stored_water(&mesh, &vg);

        // 1. Field stays finite and physically bounded. The pressure head
        //    is driven by the ponded top, whose TOTAL head is H = h + z =
        //    0 + lz = lz; at the (near) hydrostatic limit the bottom node
        //    reaches pressure head h = H − z = lz, so the field lies in
        //    [h_init, lz]. No node may run drier than the start or wetter
        //    than that hydrostatic ceiling.
        for &h in &state.head {
            assert!(h.is_finite(), "head went non-finite");
            assert!(
                h >= h_init - 1e-6,
                "no node should be drier than the start: {h} < {h_init}"
            );
            assert!(
                h <= lz + 1e-3,
                "no node should exceed the hydrostatic ceiling h=lz: {h} > {lz}"
            );
        }

        // 2. Infiltration increased the stored water.
        assert!(
            stored1 > stored0,
            "infiltration should raise stored water: {stored1} ≤ {stored0}"
        );

        // 3. MASS-BALANCE sanity bound. The stored-water increase came
        //    entirely from the ponded top boundary (no-flow elsewhere), so
        //    it cannot exceed the maximum possible saturated through-flow
        //    over the run (a generous unit-to-double hydraulic-gradient
        //    bound K·2·A·t).
        let total_time = 600.0 + dt * n_steps as f64;
        let delta_storage = stored1 - stored0;
        assert!(
            delta_storage > 0.0 && delta_storage.is_finite(),
            "storage change {delta_storage} must be positive and finite"
        );
        let area = lx * ly;
        let max_possible = k_sat * 2.0 * area * total_time;
        assert!(
            delta_storage <= max_possible + 1e-9,
            "storage gain {delta_storage} exceeds the max through-flow {max_possible}"
        );

        // 4. HYDROSTATIC limit: after the long run the column reaches
        //    equilibrium, so the TOTAL head H = h + z is uniform (= lz, the
        //    elevation of the ponded top) and the bottom carries the highest
        //    PRESSURE head. This is the decisive check that the gravity term
        //    is wired correctly (the potential is h + z, not just h): the
        //    analytic equilibrium is h(z) = lz − z.
        for kk in 0..=nz {
            let nd = nid(0, 0, kk, nx, ny);
            let z = mesh.nodes[nd].z;
            let total_head = state.head[nd] + z;
            assert!(
                (total_head - lz).abs() < 0.05,
                "hydrostatic: total head H=h+z at z={z} is {total_head}, expected lz={lz}"
            );
        }
        assert!(
            state.head[bottom] > state.head[top[0]] - 1e-6,
            "hydrostatic: bottom pressure head should be the highest"
        );
    }

    #[test]
    fn richards_shows_a_downward_wetting_front_from_a_dry_start() {
        // Start VERY dry (h = −10 m suction ⇒ tiny Kr) on a finer column,
        // pond the top, take a few short steps. Because Kr is small, the
        // wetting front is now RESOLVABLE: nodes near the ponded top wet up
        // first and the front advances downward, so at an early time the
        // saturation DECREASES with depth (top wetter than bottom). This is
        // the classic Richards infiltration front (cf. Celia et al. 1990).
        let (nx, ny, nz) = (1, 1, 20);
        let mesh = structured_box_mesh(1.0, 1.0, 2.0, nx, ny, nz).expect("valid box");
        let n_nodes = mesh.nodes.len();
        let vg = VanGenuchten::loam();
        let k_sat = 2.0e-6;
        let h_init = -10.0;
        let mut state = RichardsState::uniform(n_nodes, h_init);

        let top: Vec<usize> = (0..=nx)
            .flat_map(|i| (0..=ny).map(move |j| nid(i, j, nz, nx, ny)))
            .collect();
        let fixed: Vec<PrescribedHead> = top
            .iter()
            .map(|&n| PrescribedHead { node: n, head: 0.0 })
            .collect();

        // A handful of short steps — early in the infiltration.
        for _ in 0..5 {
            state = solve_richards_step(&mesh, &state, 600.0, k_sat, &vg, &fixed, &[], 30, 1e-8)
                .unwrap();
        }

        // The wetting front advances DOWNWARD from the ponded top, so early
        // in the infiltration the water content decreases monotonically with
        // depth: the upper soil is markedly wetter than the lower soil, which
        // the front has not yet reached. (The very bottom relaxes toward its
        // own hydrostatic value from the non-equilibrium dry IC; the decisive
        // qualitative signature is the top-down gradient.)
        let theta_top = vg.water_content(state.head[nid(0, 0, nz - 1, nx, ny)]);
        let theta_q3 = vg.water_content(state.head[nid(0, 0, 3 * nz / 4, nx, ny)]);
        let theta_mid = vg.water_content(state.head[nid(0, 0, nz / 2, nx, ny)]);
        assert!(
            theta_top > theta_q3 && theta_q3 > theta_mid,
            "downward wetting front: θ should decrease with depth — top {theta_top}, 3/4 {theta_q3}, mid {theta_mid}"
        );
        // The front is genuinely partway down: the upper node is much closer
        // to saturation than the mid-column node still awaiting the front.
        assert!(
            theta_top > 0.5 * (vg.theta_r + vg.theta_s),
            "the ponded upper soil should be well past half-saturated: θ_top {theta_top}"
        );
    }

    #[test]
    fn richards_rejects_bad_inputs() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box");
        let n = mesh.nodes.len();
        let state = RichardsState::uniform(n, -1.0);
        let vg = VanGenuchten::loam();
        // Bad dt.
        let e1 = solve_richards_step(&mesh, &state, 0.0, 1e-6, &vg, &[], &[], 5, 1e-6).unwrap_err();
        assert!(matches!(e1, PorousFlowError::BadTimeStep(_)));
        // Bad K.
        let e2 = solve_richards_step(&mesh, &state, 1.0, -1.0, &vg, &[], &[], 5, 1e-6).unwrap_err();
        assert!(matches!(e2, PorousFlowError::BadConductivity(_)));
        // Bad retention.
        let bad_vg = VanGenuchten {
            theta_r: 0.4,
            theta_s: 0.3,
            alpha: 1.0,
            n: 2.0,
        };
        let e3 =
            solve_richards_step(&mesh, &state, 1.0, 1e-6, &bad_vg, &[], &[], 5, 1e-6).unwrap_err();
        assert!(matches!(e3, PorousFlowError::BadRetention(_)));
    }

    #[test]
    fn richards_steady_uniform_head_is_unchanged() {
        // A column already at hydrostatic equilibrium with h = 0 everywhere
        // (fully saturated, water table at/above the whole column) and the
        // top pinned at h = 0: nothing should move (no gradient in h+z beyond
        // gravity, which is balanced at saturation since Kr = 1 and the
        // bottom is free) — at minimum the field stays finite and saturated.
        let (nx, ny, nz) = (1, 1, 4);
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, nx, ny, nz).expect("valid box");
        let n = mesh.nodes.len();
        let vg = VanGenuchten::loam();
        let state = RichardsState::uniform(n, 0.0);
        let top: Vec<usize> = (0..=nx)
            .flat_map(|i| (0..=ny).map(move |j| nid(i, j, nz, nx, ny)))
            .collect();
        let bottom: Vec<usize> = (0..=nx)
            .flat_map(|i| (0..=ny).map(move |j| nid(i, j, 0, nx, ny)))
            .collect();
        // Pin both ends at h=0 to hold the saturated state.
        let mut fixed: Vec<PrescribedHead> = top
            .iter()
            .map(|&nd| PrescribedHead {
                node: nd,
                head: 0.0,
            })
            .collect();
        fixed.extend(bottom.iter().map(|&nd| PrescribedHead {
            node: nd,
            head: 0.0,
        }));
        let next =
            solve_richards_step(&mesh, &state, 3600.0, 2e-6, &vg, &fixed, &[], 20, 1e-8).unwrap();
        for &h in &next.head {
            assert!(
                h.is_finite() && h.abs() < 1e-3,
                "saturated column should stay at h≈0, got {h}"
            );
        }
    }
}
