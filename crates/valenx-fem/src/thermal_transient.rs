//! Native in-process **transient** (time-dependent) thermal FEA solver.
//!
//! ## What this is
//!
//! The time-domain companion to [`crate::thermal_solver`]. Where the steady
//! solver returns the eventual equilibrium temperature, this one marches the
//! field forward in time, so you can watch a part heat up or cool down. It
//! solves the transient heat equation
//!
//! ```text
//!   ρ·cₚ·∂T/∂t = ∇·(k ∇T) + q
//! ```
//!
//! which discretises (same linear-tet basis as the steady solver) into the
//! first-order system of ODEs
//!
//! ```text
//!   C·Ṫ + K·T = Q
//! ```
//!
//! where `K` is the global **conductivity** matrix (assembled exactly as in the
//! steady solver), `Q` the heat-load vector, and `C` the **thermal capacitance
//! (heat-capacity) matrix** — the new ingredient. `C` carries the `ρ·cₚ` thermal
//! mass that makes temperature change take time.
//!
//! ## Method
//!
//! 1. **Lumped capacitance.** Each linear tet of volume `Vₑ` owns a heat
//!    capacity `ρ·cₚ·Vₑ`, split equally to its four nodes (`ρ·cₚ·Vₑ/4` each).
//!    This row-sum ("lumped") `C` is diagonal, gives the correct total heat
//!    capacity `ρ·cₚ·V`, and keeps backward-Euler steps free of the spurious
//!    over/undershoot a consistent `C` can produce.
//! 2. **Backward Euler.** The unconditionally-stable implicit step
//!
//!    ```text
//!      (C/Δt + K)·Tₙ₊₁ = (C/Δt)·Tₙ + Q
//!    ```
//!
//!    The system matrix `C/Δt + K` is constant for a fixed `Δt`, so it is
//!    assembled and Cholesky-factorised **once** and reused for every step —
//!    each step is a back-substitution.
//! 3. **Boundary conditions.** Fixed temperatures use the same large-penalty
//!    method as the steady solver; heat loads go straight into `Q`. Unlike the
//!    steady problem, a transient solve does **not** need a Dirichlet BC — the
//!    `C/Δt` term makes the system matrix positive-definite on its own.
//!
//! ## Validated
//!
//! Against closed-form results: an insulated body under a capacity-matched heat
//! input heats as the exact linear ramp `T₀ + (P/C)·t`; a single tet with three
//! nodes clamped cold relaxes as the lumped exponential `e^(−t/τ)`, `τ = C/K`;
//! and a bar driven to a fixed-temperature gradient converges to the
//! independently-computed steady field.
//!
//! ## Honest scope
//!
//! Linear, isotropic conduction on `Tet4` meshes, lumped capacitance, constant
//! loads, backward Euler — research / preliminary-design grade, the transient
//! analogue of the steady solver's "90 % case". No temperature-dependent
//! `k(T)`/`cₚ(T)`, radiation, phase change, or convection BCs yet; for those,
//! and for production transient runs, the CalculiX / Elmer adapters remain the
//! path. Within this scope the numbers are real.

use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::factorization::CscCholesky;
use nalgebra_sparse::{CooMatrix, CscMatrix};
use thiserror::Error;

use crate::thermal_solver::{
    add_diagonal, element_conductivity, tet_block, tet_volume, FixedTemperature, HeatLoad,
};
use valenx_mesh::Mesh;

/// Errors from the transient thermal solver.
#[derive(Debug, Error)]
pub enum TransientThermalError {
    /// The supplied mesh has no [`valenx_mesh::element::ElementType::Tet4`] block.
    #[error("mesh has no Tet4 element block — the transient thermal solver needs a tetrahedral volume mesh")]
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
    /// A material property or time step was non-physical (non-finite or
    /// non-positive). `name` is the offending quantity.
    #[error("invalid {name}: must be finite and positive, got {value}")]
    BadProperty {
        /// Which quantity was bad (`"conductivity"`, `"density"`, …).
        name: &'static str,
        /// The offending value.
        value: f64,
    },
    /// The initial-temperature field length did not match the node count.
    #[error("initial temperature has {got} entries but the mesh has {expected} nodes")]
    InitialFieldMismatch {
        /// Length of the supplied initial field.
        got: usize,
        /// Required length (node count).
        expected: usize,
    },
    /// A load or boundary-condition input carried a non-finite value.
    #[error("non-finite {kind} at node {node} (NaN or infinity is not a valid load/BC)")]
    InvalidLoad {
        /// 0-based node index carrying the bad value.
        node: usize,
        /// Which input was non-finite — `"heat power"` or `"prescribed temperature"`.
        kind: &'static str,
    },
    /// The Cholesky factorisation of `C/Δt + K` failed.
    #[error("linear solve failed: transient system matrix is not positive-definite")]
    SolveFailed,
}

/// Result of a transient thermal solve: the temperature field at each recorded
/// time level, starting from the initial condition at `t = 0`.
#[derive(Clone, Debug)]
pub struct TransientThermalSolution {
    /// Time of each recorded level (s); `times[0] = 0`.
    pub times: Vec<f64>,
    /// `temperature_history[step][node]` — the nodal temperature field at each
    /// time in [`times`](Self::times). The first frame is the initial field.
    pub temperature_history: Vec<Vec<f64>>,
}

impl TransientThermalSolution {
    /// The temperature field at the final recorded time.
    pub fn final_temperature(&self) -> &[f64] {
        self.temperature_history
            .last()
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// The temperature of a single node across all recorded time levels.
    pub fn node_history(&self, node: usize) -> Vec<f64> {
        self.temperature_history
            .iter()
            .map(|frame| frame.get(node).copied().unwrap_or(f64::NAN))
            .collect()
    }
}

/// The **lumped thermal capacitance** `Cᵢ = Σ ρ·cₚ·Vₑ/4` (J/K) of each node —
/// the heat capacity owned by that node, summed over the tets it belongs to.
///
/// Useful on its own (the part's thermal mass) and the diagonal `C` the
/// transient solver integrates.
///
/// # Errors
///
/// See [`TransientThermalError`] — non-physical `density`/`specific_heat`, no
/// tet block, an empty mesh, bad connectivity, or a degenerate element.
pub fn lumped_capacitance(
    mesh: &Mesh,
    density: f64,
    specific_heat: f64,
) -> Result<Vec<f64>, TransientThermalError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(TransientThermalError::EmptyMesh);
    }
    check_positive("density", density)?;
    check_positive("specific_heat", specific_heat)?;
    let tets = tet_block(mesh).ok_or(TransientThermalError::NoTetBlock)?;
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;
    let rho_cp = density * specific_heat;
    let mut c = vec![0.0_f64; n_nodes];
    for e in 0..n_elem {
        let nodes = [
            conn[4 * e] as usize,
            conn[4 * e + 1] as usize,
            conn[4 * e + 2] as usize,
            conn[4 * e + 3] as usize,
        ];
        for &nd in &nodes {
            if nd >= n_nodes {
                return Err(TransientThermalError::BadConnectivity {
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
        let vol = tet_volume(&coords).abs();
        if vol < 1.0e-18 {
            return Err(TransientThermalError::DegenerateElement(e));
        }
        let per_node = rho_cp * vol / 4.0;
        for &nd in &nodes {
            c[nd] += per_node;
        }
    }
    Ok(c)
}

/// March a heat-conduction problem forward in time with backward Euler.
///
/// `mesh` must carry a [`valenx_mesh::element::ElementType::Tet4`] block.
/// `conductivity` is `k` (W/(m·K)); `density` is `ρ` (kg/m³); `specific_heat`
/// is `cₚ` (J/(kg·K)). `fixed` / `loads` are the Dirichlet / Neumann BCs (a
/// Dirichlet BC is *optional* here). `initial_temperature` is the field at
/// `t = 0` (one entry per node). `dt` is the time step (s) and `n_steps` the
/// number of steps to take; the returned history has `n_steps + 1` frames.
///
/// # Errors
///
/// See [`TransientThermalError`].
#[allow(clippy::too_many_arguments)]
pub fn solve_transient_thermal(
    mesh: &Mesh,
    conductivity: f64,
    density: f64,
    specific_heat: f64,
    fixed: &[FixedTemperature],
    loads: &[HeatLoad],
    initial_temperature: &[f64],
    dt: f64,
    n_steps: usize,
) -> Result<TransientThermalSolution, TransientThermalError> {
    let n_nodes = mesh.nodes.len();
    if n_nodes == 0 {
        return Err(TransientThermalError::EmptyMesh);
    }
    check_positive("conductivity", conductivity)?;
    check_positive("density", density)?;
    check_positive("specific_heat", specific_heat)?;
    check_positive("dt", dt)?;
    if initial_temperature.len() != n_nodes {
        return Err(TransientThermalError::InitialFieldMismatch {
            got: initial_temperature.len(),
            expected: n_nodes,
        });
    }

    // --- global conductivity matrix K ---
    let tets = tet_block(mesh).ok_or(TransientThermalError::NoTetBlock)?;
    let conn = &tets.connectivity;
    let n_elem = conn.len() / 4;
    let mut coo = CooMatrix::<f64>::new(n_nodes, n_nodes);
    for e in 0..n_elem {
        let nodes = [
            conn[4 * e] as usize,
            conn[4 * e + 1] as usize,
            conn[4 * e + 2] as usize,
            conn[4 * e + 3] as usize,
        ];
        for &nd in &nodes {
            if nd >= n_nodes {
                return Err(TransientThermalError::BadConnectivity {
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
        let kte = element_conductivity(&coords, conductivity)
            .ok_or(TransientThermalError::DegenerateElement(e))?;
        for a in 0..4 {
            for b in 0..4 {
                let v = kte[(a, b)];
                if v != 0.0 {
                    coo.push(nodes[a], nodes[b], v);
                }
            }
        }
    }
    let k_csc = CscMatrix::from(&coo);

    // --- lumped capacitance C and the C/Δt diagonal ---
    let c = lumped_capacitance(mesh, density, specific_heat)?;
    let c_over_dt: Vec<f64> = c.iter().map(|&ci| ci / dt).collect();

    // Base transient matrix M = K + diag(C/Δt) (always positive-definite).
    let mut m = add_diagonal(&k_csc, &c_over_dt);

    // Peak diagonal of M for penalty scaling, then fold in the Dirichlet
    // penalty so the constrained rows read ≈ penalty·T = penalty·T̄.
    let mut max_diag = 0.0_f64;
    for (r, col, v) in m.triplet_iter() {
        if r == col {
            max_diag = max_diag.max(v.abs());
        }
    }
    if max_diag <= 0.0 {
        return Err(TransientThermalError::SolveFailed);
    }
    let penalty = max_diag * 1.0e8;

    // Constant part of the right-hand side: heat loads + Dirichlet penalty.
    let mut const_rhs = DVector::<f64>::zeros(n_nodes);
    for load in loads {
        if load.node >= n_nodes {
            return Err(TransientThermalError::BadConnectivity {
                elem: usize::MAX,
                node: load.node,
                n_nodes,
            });
        }
        if !load.power.is_finite() {
            return Err(TransientThermalError::InvalidLoad {
                node: load.node,
                kind: "heat power",
            });
        }
        const_rhs[load.node] += load.power;
    }
    let mut penalty_diag = vec![0.0_f64; n_nodes];
    for fb in fixed {
        if fb.node >= n_nodes {
            return Err(TransientThermalError::BadConnectivity {
                elem: usize::MAX,
                node: fb.node,
                n_nodes,
            });
        }
        if !fb.temperature.is_finite() {
            return Err(TransientThermalError::InvalidLoad {
                node: fb.node,
                kind: "prescribed temperature",
            });
        }
        penalty_diag[fb.node] += penalty;
        const_rhs[fb.node] += penalty * fb.temperature;
    }
    m = add_diagonal(&m, &penalty_diag);

    // Reject a non-finite initial field before it poisons the march.
    for (node, &t0) in initial_temperature.iter().enumerate() {
        if !t0.is_finite() {
            return Err(TransientThermalError::InvalidLoad {
                node,
                kind: "initial temperature",
            });
        }
    }

    // --- factor once, step many ---
    let chol = CscCholesky::factor(&m).map_err(|_| TransientThermalError::SolveFailed)?;

    let mut t_n = DVector::from_row_slice(initial_temperature);
    let mut times = Vec::with_capacity(n_steps + 1);
    let mut history = Vec::with_capacity(n_steps + 1);
    times.push(0.0);
    history.push(t_n.iter().copied().collect::<Vec<f64>>());

    for s in 0..n_steps {
        let mut b = const_rhs.clone();
        for i in 0..n_nodes {
            b[i] += c_over_dt[i] * t_n[i];
        }
        let sol: DMatrix<f64> = chol.solve(&b);
        t_n = sol.column(0).into_owned();
        times.push((s + 1) as f64 * dt);
        history.push(t_n.iter().copied().collect());
    }

    Ok(TransientThermalSolution {
        times,
        temperature_history: history,
    })
}

/// Validate that a property is finite and strictly positive.
fn check_positive(name: &'static str, value: f64) -> Result<(), TransientThermalError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(TransientThermalError::BadProperty { name, value });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native_solver::structured_box_mesh;
    use crate::thermal_solver::solve_steady_thermal;
    use nalgebra::Vector3;
    use valenx_mesh::element::{ElementBlock, ElementType};
    use valenx_mesh::Mesh;

    fn nid(i: usize, j: usize, k: usize, nx: usize, ny: usize) -> usize {
        i + (nx + 1) * j + (nx + 1) * (ny + 1) * k
    }

    #[test]
    fn insulated_capacity_matched_heating_is_an_exact_linear_ramp() {
        // An insulated body (no Dirichlet, no flux out) heated so that every
        // node receives Qᵢ = Cᵢ·rate must rise uniformly at exactly `rate`:
        // C·Ṫ = Q with K·(uniform) = 0 ⇒ T(t) = T₀ + rate·t, spatially flat.
        let (nx, ny, nz) = (2, 2, 2);
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, nx, ny, nz).expect("valid box");
        let n = mesh.nodes.len();
        let (rho, cp, k) = (2700.0, 900.0, 200.0); // ~aluminium
        let cap = lumped_capacitance(&mesh, rho, cp).unwrap();
        let rate = 1.5; // K/s
        let loads: Vec<HeatLoad> = (0..n)
            .map(|i| HeatLoad {
                node: i,
                power: cap[i] * rate,
            })
            .collect();
        let init = vec![300.0; n];
        let dt = 0.1;
        let n_steps = 50;
        let sol = solve_transient_thermal(&mesh, k, rho, cp, &[], &loads, &init, dt, n_steps)
            .expect("transient solve");
        let t_end = n_steps as f64 * dt;
        let expected = 300.0 + rate * t_end;
        for &temp in sol.final_temperature() {
            assert!(
                (temp - expected).abs() < 1e-6,
                "node T {temp} vs analytic ramp {expected}"
            );
        }
    }

    #[test]
    fn single_tet_relaxes_as_the_lumped_exponential() {
        // One unit tet, three nodes clamped at 0, the fourth free and started
        // hot. The free node obeys C·Ṫ + K·T = 0 ⇒ T(t) = T₀·e^(−t/τ),
        // τ = C_free / K_free,free — the textbook lumped first-order decay.
        let coords = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let mut mesh = Mesh::new("unit-tet");
        mesh.nodes = coords.to_vec();
        mesh.element_blocks.push(ElementBlock {
            element_type: ElementType::Tet4,
            connectivity: vec![0, 1, 2, 3],
        });
        let (rho, cp, k) = (1000.0, 1.0, 10.0);

        // τ from the same C and K the solver assembles (free node = index 3).
        let kte = element_conductivity(&coords, k).unwrap();
        let k_free = kte[(3, 3)];
        let c_free = lumped_capacitance(&mesh, rho, cp).unwrap()[3];
        let tau = c_free / k_free;

        let t0 = 100.0;
        let init = vec![0.0, 0.0, 0.0, t0];
        let fixed = [
            FixedTemperature {
                node: 0,
                temperature: 0.0,
            },
            FixedTemperature {
                node: 1,
                temperature: 0.0,
            },
            FixedTemperature {
                node: 2,
                temperature: 0.0,
            },
        ];
        let dt = tau * 0.005;
        let n_steps = 200; // t_end ≈ τ
        let sol = solve_transient_thermal(&mesh, k, rho, cp, &fixed, &[], &init, dt, n_steps)
            .expect("transient solve");
        let t_end = n_steps as f64 * dt;
        let got = sol.final_temperature()[3];
        let analytic = t0 * (-t_end / tau).exp();
        assert!(
            (got - analytic).abs() < 0.02 * t0,
            "free-node T {got} vs analytic e^(−t/τ) {analytic} (τ={tau})"
        );
        // And it must have actually decayed substantially (≈ e⁻¹).
        assert!(
            got < 0.5 * t0 && got > 0.2 * t0,
            "decay out of range: {got}"
        );
    }

    #[test]
    fn transient_converges_to_the_independent_steady_solution() {
        // Drive a bar to a hot/cold gradient and march long enough that the
        // transient field settles onto the field the *steady* solver returns.
        let (lx, ly, lz) = (4.0, 1.0, 1.0);
        let (nx, ny, nz) = (4, 1, 1);
        let mesh = structured_box_mesh(lx, ly, lz, nx, ny, nz).expect("valid box");
        let n = mesh.nodes.len();
        let (rho, cp, k) = (1000.0, 1.0, 50.0);
        let (t_cold, t_hot) = (300.0, 500.0);

        let mut fixed = Vec::new();
        for j in 0..=ny {
            for kk in 0..=nz {
                fixed.push(FixedTemperature {
                    node: nid(0, j, kk, nx, ny),
                    temperature: t_cold,
                });
                fixed.push(FixedTemperature {
                    node: nid(nx, j, kk, nx, ny),
                    temperature: t_hot,
                });
            }
        }
        let steady = solve_steady_thermal(&mesh, k, &fixed, &[]).unwrap();

        let init = vec![t_cold; n]; // start uniformly cold
        let dt = 0.5;
        let n_steps = 400; // backward Euler is unconditionally stable
        let sol = solve_transient_thermal(&mesh, k, rho, cp, &fixed, &[], &init, dt, n_steps)
            .expect("transient solve");

        for i in 0..n {
            let diff = (sol.final_temperature()[i] - steady.temperature[i]).abs();
            assert!(
                diff < 2.0,
                "node {i}: transient {} has not converged to steady {} (Δ={diff})",
                sol.final_temperature()[i],
                steady.temperature[i]
            );
        }
    }

    #[test]
    fn rejects_non_physical_time_step_and_properties() {
        let mesh = structured_box_mesh(1.0, 1.0, 1.0, 1, 1, 1).expect("valid box");
        let n = mesh.nodes.len();
        let init = vec![300.0; n];
        // Bad dt.
        let err = solve_transient_thermal(&mesh, 50.0, 1000.0, 900.0, &[], &[], &init, 0.0, 10)
            .unwrap_err();
        assert!(matches!(
            err,
            TransientThermalError::BadProperty { name: "dt", .. }
        ));
        // Bad density.
        let err = solve_transient_thermal(&mesh, 50.0, -1.0, 900.0, &[], &[], &init, 0.1, 10)
            .unwrap_err();
        assert!(matches!(
            err,
            TransientThermalError::BadProperty {
                name: "density",
                ..
            }
        ));
        // Initial-field length mismatch.
        let err = solve_transient_thermal(&mesh, 50.0, 1000.0, 900.0, &[], &[], &[300.0], 0.1, 10)
            .unwrap_err();
        assert!(matches!(
            err,
            TransientThermalError::InitialFieldMismatch { .. }
        ));
    }
}
