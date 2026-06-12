//! Extracellular electric field from a stimulating electrode in tissue.
//!
//! The quasi-static current-conduction problem `−∇·(σ∇φ) = I` is identical
//! in form to the steady heat equation `−∇·(k∇T) = q`, so this module
//! builds a tetrahedral tissue mesh and hands it to `valenx-fem`'s
//! [`solve_steady_thermal`] with the substitutions σ↔k (conductivity),
//! φ↔T (potential), and injected current↔heat source.
//!
//! Everything is solved in SI at the FEM boundary: node coordinates in
//! **metres**, conductivity in **S/m**, injected current in **amperes**, so
//! the solver's "temperature" output is the potential in **volts**, which we
//! return as **mV**.

use nalgebra::Vector3;
use valenx_fem::thermal_solver::{solve_steady_thermal, FixedTemperature, HeatLoad};
use valenx_mesh::element::{ElementBlock, ElementType};
use valenx_mesh::Mesh;

use crate::error::{NeuroError, Result};

/// A structured cubic block of tissue, centred at the origin, meshed with
/// linear tetrahedra (a conforming Kuhn 6-tet split of each grid cell).
pub struct TissueGrid {
    mesh: Mesh,
    n: usize,
    spacing_m: f64,
    sigma_s_m: f64,
}

/// Signed volume of the tet `(a,b,c,d)`; positive for a right-handed
/// vertex ordering.
fn signed_vol(p: &[Vector3<f64>], t: [u32; 4]) -> f64 {
    let a = p[t[0] as usize];
    (p[t[1] as usize] - a).dot(&((p[t[2] as usize] - a).cross(&(p[t[3] as usize] - a))))
}

/// Reorder a tet's vertices so its signed volume is positive (swap the last
/// two if it came out left-handed), keeping every element positively
/// oriented for the conductivity assembly.
fn orient_positive(p: &[Vector3<f64>], t: [u32; 4]) -> [u32; 4] {
    if signed_vol(p, t) < 0.0 {
        [t[0], t[1], t[3], t[2]]
    } else {
        t
    }
}

impl TissueGrid {
    /// Build a `side_mm` cube with `n` nodes per edge (`n ≥ 2`) and isotropic
    /// conductivity `sigma_s_m` (S/m), centred at the origin.
    pub fn cube(side_mm: f64, n: usize, sigma_s_m: f64) -> Self {
        assert!(n >= 2, "need at least 2 nodes per edge");
        let side_m = side_mm * 1.0e-3;
        let spacing = side_m / (n as f64 - 1.0);
        let half = side_m / 2.0;

        let mut nodes: Vec<Vector3<f64>> = Vec::with_capacity(n * n * n);
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    nodes.push(Vector3::new(
                        i as f64 * spacing - half,
                        j as f64 * spacing - half,
                        k as f64 * spacing - half,
                    ));
                }
            }
        }
        let idx = |i: usize, j: usize, k: usize| -> u32 { (i + n * j + n * n * k) as u32 };

        let mut conn: Vec<u32> = Vec::new();
        for k in 0..n - 1 {
            for j in 0..n - 1 {
                for i in 0..n - 1 {
                    let c = |di: usize, dj: usize, dk: usize| idx(i + di, j + dj, k + dk);
                    // Six tets sharing the cell's (0,0,0)-(1,1,1) diagonal —
                    // the Kuhn decomposition, which tiles space conformingly.
                    let tets = [
                        [c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(1, 1, 1)],
                        [c(0, 0, 0), c(1, 0, 0), c(1, 0, 1), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 1, 0), c(1, 1, 0), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 1, 0), c(0, 1, 1), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1)],
                        [c(0, 0, 0), c(0, 0, 1), c(0, 1, 1), c(1, 1, 1)],
                    ];
                    for t in tets {
                        conn.extend_from_slice(&orient_positive(&nodes, t));
                    }
                }
            }
        }

        let mut mesh = Mesh::new("tissue");
        mesh.nodes = nodes;
        let mut blk = ElementBlock::new(ElementType::Tet4);
        blk.connectivity = conn;
        mesh.element_blocks = vec![blk];
        mesh.recompute_stats();

        Self {
            mesh,
            n,
            spacing_m: spacing,
            sigma_s_m,
        }
    }

    fn node_index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.n * j + self.n * self.n * k
    }

    fn center_node(&self) -> usize {
        let c = self.n / 2;
        self.node_index(c, c, c)
    }

    /// Every node on an outer face, as a Dirichlet BC at `temperature`.
    fn boundary_fixed(&self, temperature: f64) -> Vec<FixedTemperature> {
        let n = self.n;
        let mut out = Vec::new();
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    if i == 0 || i == n - 1 || j == 0 || j == n - 1 || k == 0 || k == n - 1 {
                        out.push(FixedTemperature {
                            node: self.node_index(i, j, k),
                            temperature,
                        });
                    }
                }
            }
        }
        out
    }

    /// Solve the field for arbitrary Dirichlet (`fixed`, in volts) and nodal
    /// current (`loads`, in amperes) boundary conditions; returns the nodal
    /// potential in **mV**.
    fn solve_raw(&self, fixed: &[FixedTemperature], loads: &[HeatLoad]) -> Result<Vec<f64>> {
        let sol = solve_steady_thermal(&self.mesh, self.sigma_s_m, fixed, loads)
            .map_err(|e| NeuroError::Solver(e.to_string()))?;
        Ok(sol.temperature.iter().map(|v| v * 1.0e3).collect())
    }

    /// Inject `current_ua` (µA) at the central node, ground the outer faces,
    /// and solve the extracellular potential field.
    pub fn solve_point_source(&self, current_ua: f64) -> Result<ExtracellularField> {
        let loads = vec![HeatLoad {
            node: self.center_node(),
            power: current_ua * 1.0e-6,
        }];
        let fixed = self.boundary_fixed(0.0);
        let potential_mv = self.solve_raw(&fixed, &loads)?;
        Ok(ExtracellularField {
            potential_mv,
            n: self.n,
            spacing_m: self.spacing_m,
        })
    }

    /// Number of nodes per edge.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Grid node spacing (metres).
    pub fn spacing_m(&self) -> f64 {
        self.spacing_m
    }

    /// Steady **conductive** temperature rise (K, per node) from a `power_w`
    /// point heat source at the centre, with thermal conductivity `k_w_per_mk`
    /// (W/m·K) and the outer boundary held at the baseline (ΔT = 0). Reuses
    /// the same FEM solver as the electric field (`−∇·(k∇T)=q`). Pure
    /// conduction; blood-perfusion cooling is future work (see RFC 0011).
    pub fn solve_heat_raw(&self, k_w_per_mk: f64, power_w: f64) -> Result<Vec<f64>> {
        let loads = vec![HeatLoad {
            node: self.center_node(),
            power: power_w,
        }];
        let fixed = self.boundary_fixed(0.0);
        let sol = solve_steady_thermal(&self.mesh, k_w_per_mk, &fixed, &loads)
            .map_err(|e| NeuroError::Solver(e.to_string()))?;
        Ok(sol.temperature)
    }
}

/// A solved extracellular potential field over a [`TissueGrid`].
pub struct ExtracellularField {
    /// Per-node extracellular potential (mV), indexed like the grid
    /// (`i + n·j + n²·k`).
    pub potential_mv: Vec<f64>,
    n: usize,
    spacing_m: f64,
}

impl ExtracellularField {
    fn node_index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.n * j + self.n * self.n * k
    }

    /// Potential (mV) at the point `pos` (metres; the grid is centred at the
    /// origin), trilinearly interpolated from the eight surrounding nodes.
    pub fn potential_mv_at(&self, pos: Vector3<f64>) -> f64 {
        let c = (self.n as f64 - 1.0) / 2.0;
        let nmax = self.n as f64 - 1.0;
        let g = |coord: f64| (coord / self.spacing_m + c).clamp(0.0, nmax);
        let (gx, gy, gz) = (g(pos.x), g(pos.y), g(pos.z));
        let (i0, j0, k0) = (
            gx.floor() as usize,
            gy.floor() as usize,
            gz.floor() as usize,
        );
        let (i1, j1, k1) = (
            (i0 + 1).min(self.n - 1),
            (j0 + 1).min(self.n - 1),
            (k0 + 1).min(self.n - 1),
        );
        let (fx, fy, fz) = (gx - i0 as f64, gy - j0 as f64, gz - k0 as f64);
        let p = |i: usize, j: usize, k: usize| self.potential_mv[self.node_index(i, j, k)];
        let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
        let c00 = lerp(p(i0, j0, k0), p(i1, j0, k0), fx);
        let c10 = lerp(p(i0, j1, k0), p(i1, j1, k0), fx);
        let c01 = lerp(p(i0, j0, k1), p(i1, j0, k1), fx);
        let c11 = lerp(p(i0, j1, k1), p(i1, j1, k1), fx);
        let c0 = lerp(c00, c10, fy);
        let c1 = lerp(c01, c11, fy);
        lerp(c0, c1, fz)
    }

    /// Potential (mV) at distance `r_mm` along +x from the grid centre.
    pub fn potential_mv_at_radius_x(&self, r_mm: f64) -> f64 {
        self.potential_mv_at(Vector3::new(r_mm * 1.0e-3, 0.0, 0.0))
    }

    /// Grid node spacing (metres) — the natural scale at which to sample the
    /// field's curvature.
    pub fn node_spacing_m(&self) -> f64 {
        self.spacing_m
    }
}

/// Analytic potential (mV) of a point current source `current_ua` (µA) in an
/// infinite homogeneous medium of conductivity `sigma_s_m` (S/m) at distance
/// `r_mm` (mm): `φ = I / (4πσr)`.
pub fn analytic_point_source_mv(current_ua: f64, sigma_s_m: f64, r_mm: f64) -> f64 {
    let i_amp = current_ua * 1.0e-6;
    let r_m = r_mm * 1.0e-3;
    (i_amp / (4.0 * std::f64::consts::PI * sigma_s_m * r_m)) * 1.0e3
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_gradient_solves_to_precision() {
        // No source; opposite x-faces clamped to 0.1 V and 0 V. A linear
        // field is represented exactly by linear tets, so the centre plane
        // sits at 50 mV to solver precision — validates the FEM solve.
        let grid = TissueGrid::cube(20.0, 11, 0.2);
        let n = 11;
        let mut fixed = Vec::new();
        for k in 0..n {
            for j in 0..n {
                fixed.push(FixedTemperature {
                    node: grid.node_index(0, j, k),
                    temperature: 0.1,
                });
                fixed.push(FixedTemperature {
                    node: grid.node_index(n - 1, j, k),
                    temperature: 0.0,
                });
            }
        }
        let phi_mv = grid.solve_raw(&fixed, &[]).expect("solve");
        let center = grid.node_index(n / 2, n / 2, n / 2);
        assert!(
            (phi_mv[center] - 50.0).abs() < 0.5,
            "linear gradient → centre ≈ 50 mV; got {}",
            phi_mv[center]
        );
    }

    #[test]
    fn point_source_field_falls_off_with_distance() {
        let grid = TissueGrid::cube(40.0, 21, 0.2);
        let field = grid.solve_point_source(100.0).expect("solve");
        let p4 = field.potential_mv_at_radius_x(4.0);
        let p8 = field.potential_mv_at_radius_x(8.0);
        assert!(
            p4 > p8 && p8 > 0.0,
            "φ must decrease with r: p4={p4} p8={p8}"
        );
        // 1/r law (loose — coarse grid + grounded finite box distort it).
        assert!(
            (1.5..3.0).contains(&(p4 / p8)),
            "1/r-like falloff expected; φ(4)/φ(8) = {}",
            p4 / p8
        );
        // Right order of magnitude vs the infinite-medium I/(4πσr).
        let ana4 = analytic_point_source_mv(100.0, 0.2, 4.0);
        assert!(
            (0.4..2.0).contains(&(p4 / ana4)),
            "φ within mesh error of analytic; φ(4)={p4} analytic={ana4}"
        );
    }
}
