//! Tissue heating from stimulation power (Pennes bioheat; v1 = conduction).
//!
//! v1 solves steady **heat conduction** from the deposited stimulation power,
//! reusing the same FEM machinery as the extracellular field — the operators
//! `−∇·(k∇T)=q` and `−∇·(σ∇φ)=I` are identical. Blood-perfusion cooling — the
//! term that turns conduction into the full Pennes bioheat equation — is
//! future work (see RFC 0011).

use crate::error::Result;
use crate::field::TissueGrid;

/// A solved steady temperature-rise field (K) over a [`TissueGrid`].
pub struct BioheatField {
    /// Per-node temperature rise above baseline (K), grid-indexed
    /// (`i + n·j + n²·k`).
    pub delta_t_k: Vec<f64>,
    n: usize,
    spacing_m: f64,
}

impl BioheatField {
    /// Temperature rise (K) at the node nearest distance `r_mm` along +x from
    /// the grid centre.
    pub fn delta_t_k_at_radius_x(&self, r_mm: f64) -> f64 {
        let c = self.n / 2;
        let di = (r_mm * 1.0e-3 / self.spacing_m).round() as usize;
        let i = (c + di).min(self.n - 1);
        self.delta_t_k[i + self.n * c + self.n * self.n * c]
    }
}

/// Solve the steady conductive temperature rise from a `power_w`-watt point
/// heat source at the electrode, in tissue of thermal conductivity
/// `k_w_per_mk` (W/m·K).
pub fn solve_point_heat(grid: &TissueGrid, k_w_per_mk: f64, power_w: f64) -> Result<BioheatField> {
    let delta_t_k = grid.solve_heat_raw(k_w_per_mk, power_w)?;
    Ok(BioheatField {
        delta_t_k,
        n: grid.n(),
        spacing_m: grid.spacing_m(),
    })
}

/// Analytic temperature rise (K) of a point heat source `power_w` (W) in an
/// infinite homogeneous medium of thermal conductivity `k_w_per_mk` (W/m·K)
/// at radius `r_mm`: `ΔT = Q / (4πk r)`.
pub fn analytic_point_heat_k(power_w: f64, k_w_per_mk: f64, r_mm: f64) -> f64 {
    power_w / (4.0 * std::f64::consts::PI * k_w_per_mk * (r_mm * 1.0e-3))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_heat_source_matches_analytic() {
        let grid = TissueGrid::cube(40.0, 21, 0.2); // σ is unused for the thermal solve
        let k = 0.5; // brain thermal conductivity (W/m·K)
        let field = solve_point_heat(&grid, k, 0.01).expect("solve"); // 10 mW source
        let num4 = field.delta_t_k_at_radius_x(4.0);
        let num8 = field.delta_t_k_at_radius_x(8.0);
        assert!(num4 > num8 && num8 > 0.0, "ΔT must fall with distance: {num4} {num8}");
        let ana4 = analytic_point_heat_k(0.01, k, 4.0);
        assert!(
            (0.4..2.0).contains(&(num4 / ana4)),
            "ΔT within mesh error of Q/(4πk r); ΔT(4mm)={num4} analytic={ana4}"
        );
    }
}
