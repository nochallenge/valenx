//! In-house **2-D SIMP topology optimization** — generative structural design.
//!
//! Given a rectangular design domain, a load case (fixed degrees of freedom +
//! a point load) and a target material budget, the solver iteratively strips
//! material away until what is left is (close to) the **minimum-compliance**
//! (stiffest-for-the-weight) structure. The result is the characteristic
//! organic, truss-like layout generative design is known for — the same family
//! of result as Sigmund's classic "99-line" educational SIMP code (O. Sigmund,
//! *A 99 line topology optimization code written in Matlab*, Struct. Multidisc.
//! Optim. 21, 2001), reimplemented here from scratch in Rust.
//!
//! ## Method (SIMP / density-based topology optimization)
//!
//! The domain is discretised into an `nx × ny` grid of bilinear quadrilateral
//! (**Q4**) plane-stress finite elements. Each element `e` carries a *density*
//! design variable `x_e ∈ [x_min, 1]`. The **SIMP** ("Solid Isotropic Material
//! with Penalization") interpolation maps density to Young's modulus
//!
//! ```text
//! E_e(x_e) = E_min + x_e^p · (E_0 − E_min)
//! ```
//!
//! with penalty `p = 3` (the standard value) pushing intermediate densities
//! toward 0/1 so the design becomes near-discrete (solid/void). The objective
//! is the structural **compliance** `c = fᵀu = uᵀK u` (minimising compliance =
//! maximising stiffness), subject to a **volume constraint** `mean(x) ≤ vol_frac`.
//!
//! Each iteration:
//! 1. **Assemble** the global stiffness `K(x)` (sparse, SPD after the fixed
//!    DOFs are pinned) from the per-element `E_e(x_e) · k0`.
//! 2. **Solve** `K u = f` for the displacement `u` using the valenx-fem
//!    **faer** sparse-Cholesky backend ([`valenx_fem::faer_solver::solve_spd`]).
//! 3. **Element compliance & sensitivity**: `c_e = u_eᵀ k0 u_e` and
//!    `dc/dx_e = −p · x_e^(p−1) · (E_0 − E_min) · c_e` (≤ 0 — adding material
//!    can only reduce compliance).
//! 4. **Filter** the sensitivities with a linear "hat" convolution of radius
//!    `r_min` (the standard mesh-independence / checkerboard-suppression
//!    filter).
//! 5. **Optimality-Criteria (OC)** update: a Lagrange-multiplier bisection
//!    finds the material price that meets the volume target, and each density
//!    moves by the OC rule (clamped by a move limit).
//!
//! The loop stops when the largest density change between iterations falls
//! below `tol` or `max_iter` is reached. Output is the optimised density field
//! plus the per-iteration compliance history.
//!
//! ## Properties (validated in the tests)
//!
//! * The MBB-beam and cantilever benchmarks converge to the known organic
//!   truss layouts (material concentrated on the load paths, void elsewhere).
//! * The achieved volume fraction meets the target to within a few percent.
//! * Compliance decreases monotonically-ish then plateaus.
//! * The result is **deterministic** for fixed inputs (no randomness).
//!
//! The struct is framework-free (no `egui`, no app state) so it is cheap to
//! unit-test and reuse; the workbench (`valenx-app::topopt_workbench`) owns the
//! GUI, animation and heat-map rendering.

use nalgebra::DVector;
use nalgebra_sparse::{CooMatrix, CscMatrix};
use valenx_fem::faer_solver::{default_backend_for, solve_spd};

/// Young's modulus of fully-solid material (`x = 1`). The compliance scales
/// linearly with `1/E_0`, so the absolute value is just a unit choice; `1.0`
/// is the conventional normalised stiffness.
pub const E0: f64 = 1.0;
/// Young's modulus floor for void material (`x = 0`). A small positive value
/// (Ersatz material) keeps the global stiffness non-singular so the FE solve
/// never sees a free-floating void region.
pub const E_MIN: f64 = 1e-9;
/// Poisson's ratio of the (isotropic) base material. `0.3` is the textbook
/// value used by the classic 99-line benchmark.
pub const NU: f64 = 0.3;
/// SIMP move limit: the maximum change in any density per OC update. `0.2` is
/// the standard value (keeps the update stable).
const MOVE_LIMIT: f64 = 0.2;
/// OC update damping exponent (`η` in the optimality-criteria rule). `0.5` is
/// the standard value.
const OC_DAMP: f64 = 0.5;

/// Which canonical load case to optimise. Both presets use a unit point load
/// and the classic textbook supports; they differ only in geometry/symmetry,
/// which is what produces the visibly different optimal layouts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LoadCase {
    /// **MBB beam** (half-model). The Messerschmitt–Bölkow–Blohm beam: a
    /// downward point load at the **top-left** corner, a symmetry roller on the
    /// left edge (x-direction pinned) and a single vertical roller at the
    /// **bottom-right** corner. Converges to the iconic arched truss.
    MbbBeam,
    /// **Cantilever**: the whole **left edge fully clamped**, a downward point
    /// load at the **middle of the right edge**. Converges to a fanned-out
    /// load-bearing truss reaching back to the wall.
    Cantilever,
}

impl LoadCase {
    /// Short stable id (used by the workbench combo / agent bridge).
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            LoadCase::MbbBeam => "mbb",
            LoadCase::Cantilever => "cantilever",
        }
    }

    /// Human-readable label.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            LoadCase::MbbBeam => "MBB beam",
            LoadCase::Cantilever => "Cantilever",
        }
    }

    /// Parse a stable id back to a [`LoadCase`] (case-insensitive; accepts a
    /// couple of friendly aliases). `None` on an unknown id.
    #[must_use]
    pub fn from_id(s: &str) -> Option<LoadCase> {
        match s.trim().to_ascii_lowercase().as_str() {
            "mbb" | "mbbbeam" | "beam" => Some(LoadCase::MbbBeam),
            "cantilever" | "cant" => Some(LoadCase::Cantilever),
            _ => None,
        }
    }
}

/// Inputs to one SIMP run. Construct with [`TopOptParams::new`] (which clamps
/// every field to a sane range) so the optimiser cannot be handed nonsense.
#[derive(Clone, Copy, Debug)]
pub struct TopOptParams {
    /// Number of elements along the horizontal (x) axis.
    pub nx: usize,
    /// Number of elements along the vertical (y) axis.
    pub ny: usize,
    /// Target volume fraction `∈ (0, 1]` — the fraction of the domain that may
    /// remain solid (the material budget).
    pub vol_frac: f64,
    /// SIMP penalty exponent `p` (≥ 1). `3.0` is standard; higher = more
    /// black/white.
    pub penalty: f64,
    /// Sensitivity-filter radius in element units (≥ 1.0). Larger = smoother,
    /// thicker members and stronger checkerboard suppression.
    pub r_min: f64,
    /// Maximum number of optimisation iterations.
    pub max_iter: usize,
    /// Which load case to solve.
    pub load_case: LoadCase,
}

impl TopOptParams {
    /// Build a parameter set, clamping every field to a safe range:
    /// `nx, ny ∈ [4, 200]`, `vol_frac ∈ [0.05, 1.0]`, `penalty ∈ [1, 8]`,
    /// `r_min ∈ [1.0, 16.0]`, `max_iter ∈ [1, 1000]`.
    #[must_use]
    pub fn new(
        nx: usize,
        ny: usize,
        vol_frac: f64,
        penalty: f64,
        r_min: f64,
        max_iter: usize,
        load_case: LoadCase,
    ) -> Self {
        Self {
            nx: nx.clamp(4, 200),
            ny: ny.clamp(4, 200),
            vol_frac: vol_frac.clamp(0.05, 1.0),
            penalty: penalty.clamp(1.0, 8.0),
            r_min: r_min.clamp(1.0, 16.0),
            max_iter: max_iter.clamp(1, 1000),
            load_case,
        }
    }
}

impl Default for TopOptParams {
    fn default() -> Self {
        // The canonical 60×20 MBB beam at 50 % volume, p=3, r_min=1.5.
        Self::new(60, 20, 0.5, 3.0, 1.5, 60, LoadCase::MbbBeam)
    }
}

/// The result of a SIMP run: the optimised density field plus the convergence
/// history.
#[derive(Clone, Debug)]
pub struct TopOptResult {
    /// Element grid width (x).
    pub nx: usize,
    /// Element grid height (y).
    pub ny: usize,
    /// Optimised element densities, **row-major** (`idx = y·nx + x`), length
    /// `nx·ny`, each in `[x_min, 1]`. `1` = solid (black), `0` = void (white).
    pub density: Vec<f64>,
    /// Compliance `c = fᵀu` recorded **once per iteration** (the objective
    /// being minimised). Should decrease then plateau.
    pub compliance_history: Vec<f64>,
    /// Number of iterations actually run.
    pub iterations: usize,
    /// Final achieved volume fraction `mean(density)`.
    pub volume_fraction: f64,
    /// Per-iteration density-field snapshots (row-major, same layout as
    /// [`Self::density`]), captured **only** when the run requested them (see
    /// [`run_with_snapshots`]). Empty for a plain [`run`]. The workbench uses
    /// these to *animate* the structure forming iteration by iteration. The
    /// last snapshot equals [`Self::density`].
    pub snapshots: Vec<Vec<f64>>,
}

impl TopOptResult {
    /// Density at element `(x, y)` (row-major), clamped to the grid bounds.
    #[inline]
    #[must_use]
    pub fn density_at(&self, x: usize, y: usize) -> f64 {
        if self.nx == 0 || self.ny == 0 {
            return 0.0;
        }
        let x = x.min(self.nx - 1);
        let y = y.min(self.ny - 1);
        self.density[y * self.nx + x]
    }

    /// The final (last recorded) compliance, or `0.0` if no iteration ran.
    #[must_use]
    pub fn final_compliance(&self) -> f64 {
        self.compliance_history.last().copied().unwrap_or(0.0)
    }
}

// ---------------------------------------------------------------------------
// Q4 plane-stress element stiffness (unit element, unit thickness, E factored)
// ---------------------------------------------------------------------------

/// The 8×8 stiffness matrix `k0` of a **unit-square Q4 bilinear plane-stress
/// element** with `E = 1`, unit thickness and Poisson's ratio [`NU`]. Per-
/// element stiffness is then simply `E_e · k0`, so `k0` is computed once.
///
/// This is the closed-form element matrix used by Sigmund's 99-line code
/// (derived by analytic integration of the bilinear element). DOF ordering is
/// `[u1 v1 u2 v2 u3 v3 u4 v4]` for the four corner nodes in counter-clockwise
/// order starting bottom-left.
#[must_use]
pub fn element_stiffness() -> [[f64; 8]; 8] {
    let nu = NU;
    // Sigmund's 8-vector `k`, scaled by E/(1-nu^2).
    let k = [
        0.5 - nu / 6.0,
        0.125 + nu / 8.0,
        -0.25 - nu / 12.0,
        -0.125 + 3.0 * nu / 8.0,
        -0.25 + nu / 12.0,
        -0.125 - nu / 8.0,
        nu / 6.0,
        0.125 - 3.0 * nu / 8.0,
    ];
    let f = E0 / (1.0 - nu * nu);
    // The symmetric pattern of indices into `k` (Sigmund's KE assembly).
    #[rustfmt::skip]
    let idx: [[usize; 8]; 8] = [
        [0, 1, 2, 3, 4, 5, 6, 7],
        [1, 0, 7, 6, 5, 4, 3, 2],
        [2, 7, 0, 5, 6, 3, 4, 1],
        [3, 6, 5, 0, 7, 2, 1, 4],
        [4, 5, 6, 7, 0, 1, 2, 3],
        [5, 4, 3, 2, 1, 0, 7, 6],
        [6, 3, 4, 1, 2, 7, 0, 5],
        [7, 2, 1, 4, 3, 6, 5, 0],
    ];
    let mut ke = [[0.0f64; 8]; 8];
    for (i, row) in idx.iter().enumerate() {
        for (j, &m) in row.iter().enumerate() {
            ke[i][j] = f * k[m];
        }
    }
    ke
}

// ---------------------------------------------------------------------------
// DOF / node numbering
// ---------------------------------------------------------------------------

/// Node grid is `(nx+1) × (ny+1)`; node `(i, j)` (i along x, j along y) has
/// global index `j·(nx+1) + i`, and DOFs `2·idx` (x) and `2·idx + 1` (y).
#[inline]
fn node_id(nx: usize, i: usize, j: usize) -> usize {
    j * (nx + 1) + i
}

/// The eight global DOF indices of element `(ex, ey)`, in the same corner
/// order as [`element_stiffness`] (CCW from bottom-left).
#[inline]
fn element_dofs(nx: usize, ex: usize, ey: usize) -> [usize; 8] {
    let n0 = node_id(nx, ex, ey); // bottom-left
    let n1 = node_id(nx, ex + 1, ey); // bottom-right
    let n2 = node_id(nx, ex + 1, ey + 1); // top-right
    let n3 = node_id(nx, ex, ey + 1); // top-left
    [
        2 * n0,
        2 * n0 + 1,
        2 * n1,
        2 * n1 + 1,
        2 * n2,
        2 * n2 + 1,
        2 * n3,
        2 * n3 + 1,
    ]
}

// ---------------------------------------------------------------------------
// Boundary conditions / loads per case
// ---------------------------------------------------------------------------

/// The fixed (pinned) global DOFs and the load vector for a case on an
/// `nx × ny` element grid. The load is a unit downward (`−1` on the y-DOF)
/// point load at the case's load node.
fn boundary_conditions(params: &TopOptParams) -> (Vec<usize>, DVector<f64>) {
    let (nx, ny) = (params.nx, params.ny);
    let n_dof = 2 * (nx + 1) * (ny + 1);
    let mut f = DVector::zeros(n_dof);
    let mut fixed: Vec<usize> = Vec::new();

    match params.load_case {
        LoadCase::MbbBeam => {
            // Symmetry: the whole LEFT edge has its x-DOF pinned (the half-beam
            // symmetry plane).
            for j in 0..=ny {
                fixed.push(2 * node_id(nx, 0, j)); // u_x = 0 on x=0
            }
            // A single vertical roller at the bottom-right corner: y-DOF pinned.
            fixed.push(2 * node_id(nx, nx, 0) + 1);
            // Downward unit load at the top-left corner (the MBB load point).
            f[2 * node_id(nx, 0, ny) + 1] = -1.0;
        }
        LoadCase::Cantilever => {
            // The whole LEFT edge fully clamped (both DOFs).
            for j in 0..=ny {
                let n = node_id(nx, 0, j);
                fixed.push(2 * n);
                fixed.push(2 * n + 1);
            }
            // Downward unit load at the middle of the RIGHT edge.
            let jmid = ny / 2;
            f[2 * node_id(nx, nx, jmid) + 1] = -1.0;
        }
    }

    fixed.sort_unstable();
    fixed.dedup();
    (fixed, f)
}

// ---------------------------------------------------------------------------
// Sensitivity filter (linear hat convolution, radius r_min)
// ---------------------------------------------------------------------------

/// Precomputed sparse filter weights: for each element, the list of
/// `(neighbour_index, weight)` pairs within `r_min`, with `weight = r_min − dist`
/// (the standard linear "cone"/hat kernel), plus the per-element weight sums.
struct Filter {
    /// For element `e`, the neighbours `(idx, weight)` inside `r_min`.
    neigh: Vec<Vec<(usize, f64)>>,
    /// Sum of weights for each element (the convolution normaliser).
    wsum: Vec<f64>,
}

impl Filter {
    /// Build the filter for an `nx × ny` grid and radius `r_min` (element units).
    fn build(nx: usize, ny: usize, r_min: f64) -> Self {
        let n = nx * ny;
        let mut neigh = vec![Vec::new(); n];
        let mut wsum = vec![0.0f64; n];
        let r = r_min.ceil() as isize;
        for ey in 0..ny {
            for ex in 0..nx {
                let e = ey * nx + ex;
                let mut sum = 0.0;
                for dy in -r..=r {
                    for dx in -r..=r {
                        let kx = ex as isize + dx;
                        let ky = ey as isize + dy;
                        if kx < 0 || ky < 0 || kx >= nx as isize || ky >= ny as isize {
                            continue;
                        }
                        let dist = ((dx * dx + dy * dy) as f64).sqrt();
                        let w = r_min - dist;
                        if w > 0.0 {
                            let k = ky as usize * nx + kx as usize;
                            neigh[e].push((k, w));
                            sum += w;
                        }
                    }
                }
                wsum[e] = sum;
            }
        }
        Self { neigh, wsum }
    }

    /// Apply the classic **sensitivity filter** (Sigmund 1997/2001): replace
    /// each element sensitivity by a density-weighted, distance-weighted
    /// average of its neighbourhood. This makes the design mesh-independent and
    /// removes checkerboarding. `x` is the density field, `dc` the raw
    /// sensitivities; returns the filtered sensitivities.
    fn apply(&self, x: &[f64], dc: &[f64]) -> Vec<f64> {
        let n = x.len();
        let mut filtered = vec![0.0f64; n];
        for e in 0..n {
            let mut acc = 0.0;
            for &(k, w) in &self.neigh[e] {
                acc += w * x[k] * dc[k];
            }
            // Divide by (weight-sum · max(x_e, γ)) per the sensitivity filter.
            filtered[e] = acc / (self.wsum[e] * x[e].max(1e-3));
        }
        filtered
    }
}

// ---------------------------------------------------------------------------
// Global assembly + FE solve (faer backend)
// ---------------------------------------------------------------------------

/// Assemble the global stiffness `K(x)` in CSC form, with the fixed DOFs
/// enforced by the **large-penalty** (stiffened-diagonal) method so the system
/// stays SPD and the same size — exactly the boundary treatment the valenx-fem
/// faer path was built for. `ke` is the unit-`E` element matrix; `x` the
/// densities; `simp` maps density → `E` factor.
fn assemble(
    params: &TopOptParams,
    ke: &[[f64; 8]; 8],
    x: &[f64],
    fixed: &[usize],
) -> CscMatrix<f64> {
    let (nx, ny) = (params.nx, params.ny);
    let n_dof = 2 * (nx + 1) * (ny + 1);
    let p = params.penalty;
    let mut coo = CooMatrix::<f64>::new(n_dof, n_dof);
    for ey in 0..ny {
        for ex in 0..nx {
            let xe = x[ey * nx + ex];
            let e_factor = E_MIN + xe.powf(p) * (E0 - E_MIN);
            let dofs = element_dofs(nx, ex, ey);
            for (a, &ga) in dofs.iter().enumerate() {
                for (b, &gb) in dofs.iter().enumerate() {
                    let v = e_factor * ke[a][b];
                    if v != 0.0 {
                        coo.push(ga, gb, v); // CooMatrix accumulates duplicates
                    }
                }
            }
        }
    }
    let mut csc = CscMatrix::from(&coo);
    // Penalty boundary conditions: add a huge stiffness on each fixed DOF's
    // diagonal so its displacement is driven to ~0 while K stays SPD.
    let big = 1.0e12;
    for &d in fixed {
        // Find the (d, d) entry and add to it; CSC has it because every node
        // touches at least one element (diagonal is structurally present).
        add_to_diagonal(&mut csc, d, big);
    }
    csc
}

/// Add `value` to the `(d, d)` diagonal entry of a CSC matrix in place. The
/// diagonal is always structurally present here (every DOF belongs to ≥1
/// element), so this is a direct value bump, no re-symbolic-ation.
fn add_to_diagonal(csc: &mut CscMatrix<f64>, d: usize, value: f64) {
    let (col_offsets, row_indices, values) = csc.csc_data_mut();
    let start = col_offsets[d];
    let end = col_offsets[d + 1];
    for k in start..end {
        if row_indices[k] == d {
            values[k] += value;
            return;
        }
    }
}

// ---------------------------------------------------------------------------
// Optimality-Criteria update
// ---------------------------------------------------------------------------

/// One **Optimality-Criteria (OC)** density update toward the volume target.
/// A Lagrange multiplier `λ` (material price) is found by bisection so that the
/// updated mean density equals `vol_frac`; each density moves by the OC rule
///
/// ```text
/// x_new = clamp( x · (−dc / λ)^η ,  x ± move ,  [x_min, 1] )
/// ```
///
/// `dc` are the (filtered, ≤ 0) sensitivities. Returns the new field and the
/// largest per-element change (the convergence measure).
fn oc_update(x: &[f64], dc: &[f64], vol_frac: f64) -> (Vec<f64>, f64) {
    let n = x.len();
    let x_min = 1e-3;
    let mut lo = 1e-9;
    let mut hi = 1e9;
    let mut x_new = vec![0.0f64; n];
    // Bisection on the Lagrange multiplier until the volume target is hit.
    while (hi - lo) / (lo + hi) > 1e-4 {
        let lmid = 0.5 * (lo + hi);
        let mut vol = 0.0;
        for e in 0..n {
            // B_e = −dc/λ ≥ 0 (dc ≤ 0). The OC multiplicative step.
            let be = (-dc[e] / lmid).max(0.0);
            let step = x[e] * be.powf(OC_DAMP);
            let upper = (x[e] + MOVE_LIMIT).min(1.0);
            let lower = (x[e] - MOVE_LIMIT).max(x_min);
            x_new[e] = step.clamp(lower, upper);
            vol += x_new[e];
        }
        // Too much material at this price → raise the price (raise λ).
        if vol > vol_frac * n as f64 {
            lo = lmid;
        } else {
            hi = lmid;
        }
    }
    let change = x
        .iter()
        .zip(&x_new)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    (x_new, change)
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// Run the full SIMP loop and return the optimised density field plus the
/// compliance history. Deterministic for fixed `params` (no randomness): the
/// field is initialised uniformly at `vol_frac` and every step is a
/// closed-form computation.
///
/// Convergence tolerance is fixed at `0.01` on the max density change (the
/// 99-line default); the loop also stops at `params.max_iter`. The returned
/// [`TopOptResult::snapshots`] is empty; use [`run_with_snapshots`] to capture
/// the per-iteration fields for animation.
#[must_use]
pub fn run(params: &TopOptParams) -> TopOptResult {
    run_core(params, false)
}

/// Like [`run`], but also records the density field **after every iteration**
/// into [`TopOptResult::snapshots`] (the last snapshot equals the final
/// `density`). The workbench plays these back to animate the structure
/// emerging. Identical optimisation to [`run`] — only the bookkeeping differs.
#[must_use]
pub fn run_with_snapshots(params: &TopOptParams) -> TopOptResult {
    run_core(params, true)
}

/// Shared SIMP driver behind [`run`] / [`run_with_snapshots`]; `capture`
/// toggles the per-iteration snapshot recording.
#[must_use]
fn run_core(params: &TopOptParams, capture: bool) -> TopOptResult {
    let (nx, ny) = (params.nx, params.ny);
    let n_el = nx * ny;
    let ke = element_stiffness();
    let (fixed, f) = boundary_conditions(params);
    let filter = Filter::build(nx, ny, params.r_min);
    let n_dof = f.len();
    let backend = default_backend_for(n_dof);

    // Uniform initial density = the target volume fraction (a feasible start).
    let mut x = vec![params.vol_frac; n_el];
    let mut compliance_history: Vec<f64> = Vec::new();
    let mut snapshots: Vec<Vec<f64>> = Vec::new();
    let tol = 0.01;
    let mut iterations = 0;

    for _ in 0..params.max_iter {
        iterations += 1;

        // 1. Assemble K(x), 2. solve K u = f via the faer sparse backend.
        let k = assemble(params, &ke, &x, &fixed);
        let u = match solve_spd(&k, &f, backend) {
            Ok(u) => u,
            // A singular/indefinite system (degenerate input) ends the loop
            // gracefully rather than panicking; keep whatever history we have.
            Err(_) => break,
        };

        // 3. Element compliance c_e = u_eᵀ k0 u_e and the SIMP sensitivity
        //    dc/dx_e = −p x^(p−1)(E0−Emin) c_e.  Also accumulate total
        //    compliance c = fᵀu = Σ E_e c_e.
        let p = params.penalty;
        let mut dc = vec![0.0f64; n_el];
        let mut compliance = 0.0;
        for ey in 0..ny {
            for ex in 0..nx {
                let e = ey * nx + ex;
                let dofs = element_dofs(nx, ex, ey);
                // ce = uᵀ k0 u for this element (unit-E element energy).
                let mut ce = 0.0;
                for (a, &ga) in dofs.iter().enumerate() {
                    let mut row = 0.0;
                    for (b, &gb) in dofs.iter().enumerate() {
                        row += ke[a][b] * u[gb];
                    }
                    ce += u[ga] * row;
                }
                let xe = x[e];
                let e_factor = E_MIN + xe.powf(p) * (E0 - E_MIN);
                compliance += e_factor * ce;
                dc[e] = -p * xe.powf(p - 1.0) * (E0 - E_MIN) * ce;
            }
        }
        compliance_history.push(compliance);

        // 4. Filter the sensitivities (mesh-independence / no checkerboard).
        let dc_f = filter.apply(&x, &dc);

        // 5. OC update toward the volume target.
        let (x_new, change) = oc_update(&x, &dc_f, params.vol_frac);
        x = x_new;

        if capture {
            snapshots.push(x.clone());
        }

        if change < tol {
            break;
        }
    }

    let volume_fraction = if n_el == 0 {
        0.0
    } else {
        x.iter().sum::<f64>() / n_el as f64
    };

    TopOptResult {
        nx,
        ny,
        density: x,
        compliance_history,
        iterations,
        volume_fraction,
        snapshots,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Fraction of the design occupied by "solid-ish" elements (x > 0.5).
    fn solid_fraction(r: &TopOptResult) -> f64 {
        let solid = r.density.iter().filter(|&&v| v > 0.5).count();
        solid as f64 / r.density.len() as f64
    }

    /// Mean density of the elements in a rectangular window [x0,x1) × [y0,y1).
    fn window_mean(r: &TopOptResult, x0: usize, x1: usize, y0: usize, y1: usize) -> f64 {
        let mut s = 0.0;
        let mut n = 0;
        for y in y0..y1 {
            for x in x0..x1 {
                s += r.density_at(x, y);
                n += 1;
            }
        }
        s / n as f64
    }

    #[test]
    fn element_stiffness_is_symmetric_psd_ish() {
        let ke = element_stiffness();
        // Symmetric.
        for i in 0..8 {
            for j in 0..8 {
                assert!((ke[i][j] - ke[j][i]).abs() < 1e-12, "ke not symmetric");
            }
        }
        // A rigid-body translation in x is a zero-energy mode: K · [1,0,1,0,...]
        // = 0 (rows sum to zero over the x-DOFs).
        for i in 0..8 {
            let mut row_x = 0.0;
            for j in (0..8).step_by(2) {
                row_x += ke[i][j];
            }
            assert!(row_x.abs() < 1e-9, "x rigid-body mode not in null space");
        }
    }

    #[test]
    fn params_are_clamped() {
        let p = TopOptParams::new(1, 9999, 5.0, 99.0, 0.1, 0, LoadCase::MbbBeam);
        assert_eq!(p.nx, 4);
        assert_eq!(p.ny, 200);
        assert_eq!(p.vol_frac, 1.0);
        assert_eq!(p.penalty, 8.0);
        assert!((p.r_min - 1.0).abs() < 1e-12);
        assert_eq!(p.max_iter, 1);
    }

    #[test]
    fn load_case_id_roundtrip() {
        for lc in [LoadCase::MbbBeam, LoadCase::Cantilever] {
            assert_eq!(LoadCase::from_id(lc.id()), Some(lc));
        }
        assert_eq!(LoadCase::from_id("BEAM"), Some(LoadCase::MbbBeam));
        assert_eq!(LoadCase::from_id("cant"), Some(LoadCase::Cantilever));
        assert_eq!(LoadCase::from_id("nope"), None);
    }

    #[test]
    fn volume_target_is_met() {
        // The achieved volume fraction must match the target to a few percent.
        let p = TopOptParams::new(48, 24, 0.4, 3.0, 1.5, 40, LoadCase::MbbBeam);
        let r = run(&p);
        assert!(
            (r.volume_fraction - 0.4).abs() < 0.03,
            "volume fraction {} should be ~0.4",
            r.volume_fraction
        );
    }

    #[test]
    fn compliance_decreases_then_plateaus() {
        // The objective must drop substantially from the uniform start and end
        // up lower than it began (stiffer structure for the same material).
        let p = TopOptParams::new(60, 20, 0.5, 3.0, 1.5, 60, LoadCase::MbbBeam);
        let r = run(&p);
        let h = &r.compliance_history;
        assert!(h.len() >= 5, "expected several iterations, got {}", h.len());
        let first = h[0];
        let last = *h.last().unwrap();
        assert!(
            last < first,
            "final compliance {last} should be below the initial {first}"
        );
        // Plateau: the last few iterations change little relative to the drop.
        let n = h.len();
        let tail_change = (h[n - 1] - h[n - 2]).abs();
        let total_drop = (first - last).abs();
        assert!(
            tail_change < 0.1 * total_drop + 1e-9,
            "tail step {tail_change} should be small vs total drop {total_drop}"
        );
        // All compliances finite and positive.
        assert!(h.iter().all(|c| c.is_finite() && *c > 0.0));
    }

    #[test]
    fn mbb_converges_to_a_truss_not_a_solid_block() {
        // A converged SIMP result is mostly void with material on the load
        // paths — the solid fraction should be well below 100 % and roughly
        // around the volume target, and the field should be near-discrete.
        let p = TopOptParams::new(60, 20, 0.5, 3.0, 1.5, 80, LoadCase::MbbBeam);
        let r = run(&p);
        let sf = solid_fraction(&r);
        assert!(
            (0.2..0.75).contains(&sf),
            "solid fraction {sf} should look like a truss, not a block/empty"
        );
        // Near-discreteness: few elements stuck at grey (0.2..0.8).
        let grey = r
            .density
            .iter()
            .filter(|&&v| (0.2..0.8).contains(&v))
            .count() as f64
            / r.density.len() as f64;
        assert!(
            grey < 0.5,
            "too many grey elements ({grey}); SIMP not biting"
        );
    }

    #[test]
    fn mbb_material_concentrates_along_the_top_and_bottom_chords() {
        // The MBB optimum is an arch/truss: the loaded TOP-LEFT region and the
        // supported BOTTOM region carry far more material than the mid-height
        // interior near the free (right) side.
        let p = TopOptParams::new(80, 24, 0.5, 3.0, 2.0, 90, LoadCase::MbbBeam);
        let r = run(&p);
        // Top-left load corner band vs. a mid-height interior window: the load
        // path must be denser than the interior void.
        let top_left = window_mean(&r, 0, 12, r.ny - 4, r.ny);
        let mid_interior = window_mean(&r, r.nx / 2, r.nx / 2 + 12, r.ny / 2 - 2, r.ny / 2 + 2);
        assert!(
            top_left > mid_interior,
            "loaded top-left chord ({top_left}) should be denser than mid interior ({mid_interior})"
        );
    }

    #[test]
    fn cantilever_forms_flanges_at_the_clamped_wall() {
        // For the cantilever the bending moment is largest at the clamped (left)
        // wall, so the optimum grows an I-beam-like pair of flanges there: the
        // TOP and BOTTOM bands at the wall go (near-)solid while the MID-HEIGHT
        // (neutral axis) is hollowed to void — exactly the structural truss a
        // generative design should discover.
        let p = TopOptParams::new(64, 40, 0.4, 3.0, 2.0, 80, LoadCase::Cantilever);
        let r = run(&p);
        let wall_top = window_mean(&r, 0, 8, r.ny - 6, r.ny);
        let wall_bot = window_mean(&r, 0, 8, 0, 6);
        let wall_mid = window_mean(&r, 0, 8, r.ny / 2 - 3, r.ny / 2 + 3);
        assert!(
            wall_top > 0.7 && wall_bot > 0.7,
            "wall flanges should be near-solid (top {wall_top}, bottom {wall_bot})"
        );
        assert!(
            wall_top > wall_mid && wall_bot > wall_mid,
            "wall flanges ({wall_top}/{wall_bot}) should bracket a void neutral axis ({wall_mid})"
        );
        assert!(
            wall_mid < 0.3,
            "neutral-axis mid should hollow out ({wall_mid})"
        );
        assert!((r.volume_fraction - 0.4).abs() < 0.03);
    }

    #[test]
    fn result_is_deterministic() {
        // Identical params → bit-identical density fields and history.
        let p = TopOptParams::new(40, 20, 0.5, 3.0, 1.5, 30, LoadCase::MbbBeam);
        let a = run(&p);
        let b = run(&p);
        assert_eq!(a.density, b.density, "density field must be deterministic");
        assert_eq!(
            a.compliance_history, b.compliance_history,
            "compliance history must be deterministic"
        );
    }

    #[test]
    fn densities_stay_in_bounds_and_finite() {
        let p = TopOptParams::new(50, 25, 0.45, 3.0, 1.6, 40, LoadCase::Cantilever);
        let r = run(&p);
        assert_eq!(r.density.len(), 50 * 25);
        for &d in &r.density {
            assert!(d.is_finite(), "density must be finite");
            assert!(
                (1e-3 - 1e-6..=1.0 + 1e-6).contains(&d),
                "density {d} out of bounds"
            );
        }
    }

    #[test]
    fn snapshots_capture_each_iteration_and_end_at_the_final_field() {
        // run_with_snapshots records one density field per iteration; the last
        // snapshot must equal the returned final density (so the animation ends
        // on the optimised structure), and plain run() captures none.
        let p = TopOptParams::new(40, 20, 0.5, 3.0, 1.5, 25, LoadCase::MbbBeam);
        let r = run_with_snapshots(&p);
        assert_eq!(
            r.snapshots.len(),
            r.iterations,
            "one snapshot per iteration"
        );
        assert_eq!(
            r.snapshots.last().unwrap(),
            &r.density,
            "final snapshot must equal the final density field"
        );
        assert!(r.snapshots.iter().all(|s| s.len() == 40 * 20));
        let plain = run(&p);
        assert!(
            plain.snapshots.is_empty(),
            "plain run captures no snapshots"
        );
        // Same optimisation either way.
        assert_eq!(plain.density, r.density);
    }

    #[test]
    fn uses_faer_backend_for_the_default_grid() {
        // The default 60×20 grid has 2·61·21 = 2562 DOFs, above the faer
        // threshold — confirm the solver picks faer (i.e. the result is sane and
        // the DOF count is in faer territory), guarding the "reuses the faer
        // sparse solver" contract.
        let p = TopOptParams::default();
        let n_dof = 2 * (p.nx + 1) * (p.ny + 1);
        assert!(
            n_dof >= valenx_fem::faer_solver::FAER_DOF_THRESHOLD,
            "default grid should exercise the faer backend ({n_dof} dofs)"
        );
        let r = run(&p);
        assert!(r.iterations >= 1 && r.final_compliance() > 0.0);
    }
}
