//! The activating function — how the extracellular field drives the membrane.
//!
//! Following Rattay (1986), the extracellular potential `Vₑ` sampled along a
//! fiber drives the cable through the **second spatial derivative** of `Vₑ`
//! along the fiber. The sign is the key physics: under a **cathodic**
//! electrode `Vₑ` has a minimum (a valley), so `d²Vₑ/dx² > 0` there and the
//! membrane **depolarizes** (excitation); under an **anodic** electrode `Vₑ`
//! peaks, `d²Vₑ/dx² < 0`, and the membrane **hyperpolarizes**.

use nalgebra::Vector3;

use crate::field::ExtracellularField;

/// Activating function along a fiber parallel to the x-axis at perpendicular
/// offset `y_off_mm` from the electrode: sample `Vₑ` at `n_samples` points
/// spaced `dx_mm` apart and centred on x = 0, then return the discrete second
/// difference `(Vₑ[k-1] − 2Vₑ[k] + Vₑ[k+1]) / Δx²` (mV/m²) at each interior
/// point.
pub fn activating_along_x(
    field: &ExtracellularField,
    y_off_mm: f64,
    n_samples: usize,
    dx_mm: f64,
) -> Vec<f64> {
    let dx_m = dx_mm * 1.0e-3;
    let y = y_off_mm * 1.0e-3;
    let half = (n_samples as f64 - 1.0) / 2.0;
    let ve: Vec<f64> = (0..n_samples)
        .map(|k| {
            let x = (k as f64 - half) * dx_m;
            field.potential_mv_at(Vector3::new(x, y, 0.0))
        })
        .collect();
    (1..n_samples - 1)
        .map(|k| (ve[k - 1] - 2.0 * ve[k] + ve[k + 1]) / (dx_m * dx_m))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::TissueGrid;

    #[test]
    fn cathodic_depolarizes_under_electrode_anodic_hyperpolarizes() {
        let grid = TissueGrid::cube(40.0, 21, 0.2);
        // A fiber 2 mm from the electrode, sampled along x at the 2 mm pitch.
        let cathodic = grid.solve_point_source(-100.0).expect("solve"); // current sink
        let anodic = grid.solve_point_source(100.0).expect("solve"); // current source
        let f_cath = activating_along_x(&cathodic, 2.0, 21, 2.0);
        let f_anod = activating_along_x(&anodic, 2.0, 21, 2.0);
        let mid = f_cath.len() / 2;
        assert!(
            f_cath[mid] > 0.0,
            "cathodic → depolarizing under the electrode; f={}",
            f_cath[mid]
        );
        assert!(
            f_anod[mid] < 0.0,
            "anodic → hyperpolarizing under the electrode; f={}",
            f_anod[mid]
        );
    }
}
