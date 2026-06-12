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

/// The **point-source isotherm radius** `r = Q / (4πk·ΔT)` (mm) — the radius at which the
/// steady conductive temperature rise of a `power_w`-W (Q) point heat source falls to
/// `delta_t_k` `ΔT` (K), in tissue of thermal conductivity `k_w_per_mk` (W/m·K). It is the
/// inverse of [`analytic_point_heat_k`] solved for the radius: the thermal-reach (or
/// lesion-boundary) radius of a damaging isotherm, the sizing figure for a thermal safety
/// margin around a stimulating electrode. Inputs are positive; `ΔT → 0` sends the radius
/// to infinity (the rise reaches everywhere).
pub fn point_heat_isotherm_radius_mm(power_w: f64, k_w_per_mk: f64, delta_t_k: f64) -> f64 {
    1.0e3 * power_w / (4.0 * std::f64::consts::PI * k_w_per_mk * delta_t_k)
}

/// The **point-source power for a target temperature rise** `Q = ΔT·4πk·r` (W) — the
/// deposited power that makes a point heat source's steady conductive temperature rise
/// reach `delta_t_k` `ΔT` (K) at radius `r_mm` (mm), in tissue of thermal conductivity
/// `k_w_per_mk` (W/m·K). It is the power solve of `ΔT = Q/(4πk·r)` — the third member of
/// the point-source inverse trio with [`analytic_point_heat_k`] (the rise) and
/// [`point_heat_isotherm_radius_mm`] (the radius) — the thermal-dose planning quantity
/// (how much power reaches a target rise at a given distance). `ΔT = 0` needs no power.
pub fn point_heat_power_for_delta_t_w(delta_t_k: f64, k_w_per_mk: f64, r_mm: f64) -> f64 {
    delta_t_k * 4.0 * std::f64::consts::PI * k_w_per_mk * (r_mm * 1.0e-3)
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
        assert!(
            num4 > num8 && num8 > 0.0,
            "ΔT must fall with distance: {num4} {num8}"
        );
        let ana4 = analytic_point_heat_k(0.01, k, 4.0);
        assert!(
            (0.4..2.0).contains(&(num4 / ana4)),
            "ΔT within mesh error of Q/(4πk r); ΔT(4mm)={num4} analytic={ana4}"
        );
    }

    #[test]
    fn point_heat_isotherm_radius_inverts_the_analytic_rise() {
        // (a) ROUND-TRIP threading analytic_point_heat_k (non-tautological): the radius at
        // which ΔT occurs, fed back, reproduces ΔT — and a known radius recovers itself.
        for &(p, k, dt) in &[
            (0.01_f64, 0.5_f64, 0.5_f64),
            (0.02, 0.6, 0.2),
            (0.005, 0.49, 1.0),
        ] {
            let r = point_heat_isotherm_radius_mm(p, k, dt);
            assert!(
                (analytic_point_heat_k(p, k, r) - dt).abs() <= 1e-9 * dt,
                "ΔT(r(ΔT)) = ΔT"
            );
        }
        for &(p, k, r0) in &[(0.01_f64, 0.5_f64, 4.0_f64), (0.02, 0.6, 8.0)] {
            assert!(
                (point_heat_isotherm_radius_mm(p, k, analytic_point_heat_k(p, k, r0)) - r0).abs()
                    <= 1e-9 * r0,
                "r(ΔT(r)) = r"
            );
        }

        // (b) WORKED: P = 10 mW, k = 0.5, ΔT = 0.5 K → r = 1000·P/(4πk·ΔT) = 10/π ≈ 3.183 mm.
        let r = point_heat_isotherm_radius_mm(0.01, 0.5, 0.5);
        assert!(
            (r - 10.0 / std::f64::consts::PI).abs() <= 1e-9 * r,
            "r = 10/π ≈ 3.183 mm"
        );

        // (c) MONOTONICITY: a hotter contour sits closer to the source; more power or lower
        // conductivity pushes a given contour farther out.
        assert!(
            point_heat_isotherm_radius_mm(0.01, 0.5, 1.0)
                < point_heat_isotherm_radius_mm(0.01, 0.5, 0.5),
            "hotter contour → smaller radius"
        );
        assert!(
            point_heat_isotherm_radius_mm(0.02, 0.5, 0.5)
                > point_heat_isotherm_radius_mm(0.01, 0.5, 0.5),
            "more power → larger radius"
        );
        assert!(
            point_heat_isotherm_radius_mm(0.01, 0.25, 0.5)
                > point_heat_isotherm_radius_mm(0.01, 0.5, 0.5),
            "lower conductivity → larger radius"
        );
    }

    #[test]
    fn point_heat_power_for_delta_t_completes_the_inverse_trio() {
        // (a) ROUND-TRIP threading analytic_point_heat_k (non-tautological, both directions).
        for &(dt, k, r) in &[
            (0.5_f64, 0.5_f64, 4.0_f64),
            (0.2, 0.6, 8.0),
            (1.0, 0.49, 2.0),
        ] {
            let p = point_heat_power_for_delta_t_w(dt, k, r);
            assert!(
                (analytic_point_heat_k(p, k, r) - dt).abs() <= 1e-9 * dt,
                "ΔT(Q(ΔT)) = ΔT"
            );
        }
        for &(p, k, r) in &[(0.01_f64, 0.5_f64, 4.0_f64), (0.02, 0.6, 8.0)] {
            assert!(
                (point_heat_power_for_delta_t_w(analytic_point_heat_k(p, k, r), k, r) - p).abs()
                    <= 1e-9 * p,
                "Q(ΔT(Q)) = Q"
            );
        }

        // (b) THREAD point_heat_isotherm_radius_mm (#401) (non-tautological): the power to
        // reach ΔT at the isotherm radius of a source P0 recovers P0.
        let (p0, k, dt) = (0.01_f64, 0.5_f64, 0.5_f64);
        assert!(
            (point_heat_power_for_delta_t_w(dt, k, point_heat_isotherm_radius_mm(p0, k, dt)) - p0)
                .abs()
                <= 1e-9 * p0,
            "Q at the ΔT-isotherm radius = P0"
        );

        // (c) WORKED: ΔT = 1 K, k = 0.5, r = 1 mm → Q = 4π·0.5·1·1e-3 = 2π·1e-3 ≈ 6.283 mW.
        let q = point_heat_power_for_delta_t_w(1.0, 0.5, 1.0);
        let expected = 2.0 * std::f64::consts::PI * 1.0e-3;
        assert!((q - expected).abs() <= 1e-9 * expected, "Q = 2π·1e-3 W");

        // (d) PROPORTIONALITY: Q is linear in ΔT, k, and r.
        let base = point_heat_power_for_delta_t_w(1.0, 0.5, 1.0);
        assert!(
            (point_heat_power_for_delta_t_w(2.0, 0.5, 1.0) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "∝ ΔT"
        );
        assert!(
            (point_heat_power_for_delta_t_w(1.0, 1.0, 1.0) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "∝ k"
        );
        assert!(
            (point_heat_power_for_delta_t_w(1.0, 0.5, 2.0) - 2.0 * base).abs() <= 1e-9 * 2.0 * base,
            "∝ r"
        );

        // (e) ZERO: no temperature rise needs no power.
        assert_eq!(point_heat_power_for_delta_t_w(0.0, 0.5, 1.0), 0.0);
    }
}
