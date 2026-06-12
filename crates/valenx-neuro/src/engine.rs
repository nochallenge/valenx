//! Coupled extracellular stimulation: an electrode field driving HH axons.
//!
//! The orchestrator solves the [`crate::field`] once for the electrode
//! current, samples the extracellular potential along each axon, turns it
//! into the [`crate::cable`] activating drive, integrates the Hodgkin–Huxley
//! membrane, and reports which axons fired (recruitment).

use nalgebra::Vector3;

use crate::cable::HhCable;
use crate::error::Result;
use crate::field::{ExtracellularField, TissueGrid};

/// A straight unmyelinated axon, parallel to the x-axis at perpendicular
/// distance `depth_mm` from the electrode (which sits at the grid centre).
pub struct Axon {
    /// Perpendicular distance from the electrode (mm).
    pub depth_mm: f64,
    /// Number of compartments.
    pub n_comp: usize,
    /// Compartment length (µm).
    pub dx_um: f64,
    /// Fiber radius (µm).
    pub a_um: f64,
    /// Intracellular resistivity (Ω·cm).
    pub ri_ohm_cm: f64,
}

impl Axon {
    /// A squid-like fiber (radius 238 µm, Rᵢ 35.4 Ω·cm, 101 × 100-µm
    /// compartments ≈ 1 cm long) at `depth_mm` from the electrode.
    pub fn squid_at(depth_mm: f64) -> Self {
        // 1-mm compartments: coarse enough that the explicit RK4 integrator is
        // stable in the sub-threshold (linear) regime, while still resolving
        // the ~6 mm space constant. A finer cable would need an implicit axial
        // solver (see RFC 0011, "Unresolved questions").
        Self {
            depth_mm,
            n_comp: 21,
            dx_um: 1000.0,
            a_um: 238.0,
            ri_ohm_cm: 35.4,
        }
    }
}

/// A stimulation scene: a tissue block with a central point electrode and a
/// set of axons.
pub struct Scene {
    /// The tissue grid (the electrode sits at its centre).
    pub grid: TissueGrid,
    /// The axons in the tissue.
    pub axons: Vec<Axon>,
}

impl Scene {
    /// One squid axon at `depth_mm` in a default 40 mm gray-matter block
    /// (σ = 0.2 S/m).
    pub fn single_axon(depth_mm: f64) -> Self {
        Self {
            grid: TissueGrid::cube(40.0, 21, 0.2),
            axons: vec![Axon::squid_at(depth_mm)],
        }
    }

    /// `n` squid axons spread linearly in depth over `[0.5, 0.5 + spread_mm]`.
    pub fn bundle(n: usize, spread_mm: f64) -> Self {
        let axons = (0..n)
            .map(|i| {
                let frac = if n > 1 {
                    i as f64 / (n as f64 - 1.0)
                } else {
                    0.0
                };
                Axon::squid_at(0.5 + frac * spread_mm)
            })
            .collect();
        Self {
            grid: TissueGrid::cube(40.0, 21, 0.2),
            axons,
        }
    }
}

/// Which axons fired under a given stimulus.
pub struct Recruitment {
    fired: Vec<bool>,
}

impl Recruitment {
    /// Per-axon fired flags (in scene-axon order).
    pub fn fired(&self) -> &[bool] {
        &self.fired
    }

    /// Did any axon fire?
    pub fn any_fired(&self) -> bool {
        self.fired.iter().any(|&f| f)
    }

    /// Fraction of axons recruited (0 if there are none).
    pub fn recruited_fraction(&self) -> f64 {
        if self.fired.is_empty() {
            0.0
        } else {
            self.fired.iter().filter(|&&f| f).count() as f64 / self.fired.len() as f64
        }
    }
}

/// Per-compartment extracellular activating drive (µA/cm²) for `ax`: the
/// Rattay term `1e3·(a/2Rᵢ)·∂²Vₑ/∂x²`. The second derivative uses a stencil at
/// the **field's own resolution**, so the coarse grid's piecewise-linear kinks
/// don't blow the drive up at the much finer compartment scale.
fn activating_drive_along_axon(field: &ExtracellularField, ax: &Axon) -> Vec<f64> {
    let delta_m = field.node_spacing_m();
    let delta_cm = delta_m * 100.0;
    let a_cm = ax.a_um * 1.0e-4;
    let coeff = 1.0e3 * a_cm / (2.0 * ax.ri_ohm_cm); // µA/cm² per (mV/cm²)
    let y_m = ax.depth_mm * 1.0e-3;
    let dx_m = ax.dx_um * 1.0e-6;
    let half = (ax.n_comp as f64 - 1.0) / 2.0;
    (0..ax.n_comp)
        .map(|k| {
            let x = (k as f64 - half) * dx_m;
            let vem = field.potential_mv_at(Vector3::new(x - delta_m, y_m, 0.0));
            let ve0 = field.potential_mv_at(Vector3::new(x, y_m, 0.0));
            let vep = field.potential_mv_at(Vector3::new(x + delta_m, y_m, 0.0));
            let ve_pp = (vem - 2.0 * ve0 + vep) / (delta_cm * delta_cm); // mV/cm²
            coeff * ve_pp
        })
        .collect()
}

/// Solve the field for `electrode_ua` (negative = cathodic), drive each axon
/// through the activating function, and report which fire under a 0.5 ms
/// pulse.
pub fn stimulate(scene: &Scene, electrode_ua: f64) -> Result<Recruitment> {
    let field = scene.grid.solve_point_source(electrode_ua)?;
    let mut fired = Vec::with_capacity(scene.axons.len());
    for ax in &scene.axons {
        let mut cable = HhCable::uniform(ax.n_comp, ax.dx_um, ax.a_um, ax.ri_ohm_cm);
        let drive = activating_drive_along_axon(&field, ax);
        let run = cable.stimulate_extracellular(&drive, 1.0, 0.5, 12.0, 0.01);
        let any = (0..ax.n_comp).any(|k| run.peak_time_ms(k).is_some());
        fired.push(any);
    }
    Ok(Recruitment { fired })
}

/// Recruited fraction vs cathodic current magnitude (µA). Efficient: the field
/// is linear in current, so it is solved once and each axon's activating drive
/// is scaled across the sweep. Returns `(magnitude_µA, fraction)` pairs.
pub fn recruitment_curve(scene: &Scene, current_mags_ua: &[f64]) -> Result<Vec<(f64, f64)>> {
    let field = scene.grid.solve_point_source(-1.0)?;
    let unit_drives: Vec<Vec<f64>> = scene
        .axons
        .iter()
        .map(|ax| activating_drive_along_axon(&field, ax))
        .collect();
    let n = scene.axons.len();
    let mut out = Vec::with_capacity(current_mags_ua.len());
    for &mag in current_mags_ua {
        let mut fired = 0usize;
        for (ax, unit) in scene.axons.iter().zip(&unit_drives) {
            let drive: Vec<f64> = unit.iter().map(|d| d * mag).collect();
            let mut cable = HhCable::uniform(ax.n_comp, ax.dx_um, ax.a_um, ax.ri_ohm_cm);
            let run = cable.stimulate_extracellular(&drive, 1.0, 0.5, 12.0, 0.01);
            if (0..ax.n_comp).any(|k| run.peak_time_ms(k).is_some()) {
                fired += 1;
            }
        }
        let frac = if n == 0 { 0.0 } else { fired as f64 / n as f64 };
        out.push((mag, frac));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Smallest cathodic current magnitude (µA) that recruits a fiber at
    /// `depth_mm`. The field is linear in current, so it is solved once (at a
    /// reference current) and the bisection just scales the resulting drive.
    fn threshold_ua(depth_mm: f64) -> f64 {
        let scene = Scene::single_axon(depth_mm);
        let ax = &scene.axons[0];
        let field = scene.grid.solve_point_source(-1.0).unwrap();
        let unit_drive = activating_drive_along_axon(&field, ax);
        let fires = |scale: f64| {
            let drive: Vec<f64> = unit_drive.iter().map(|d| d * scale).collect();
            let mut cable = HhCable::uniform(ax.n_comp, ax.dx_um, ax.a_um, ax.ri_ohm_cm);
            let run = cable.stimulate_extracellular(&drive, 1.0, 0.5, 12.0, 0.01);
            (0..ax.n_comp).any(|k| run.peak_time_ms(k).is_some())
        };
        let (mut lo, mut hi) = (0.0_f64, 5000.0_f64);
        for _ in 0..18 {
            let mid = 0.5 * (lo + hi);
            if fires(mid) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }

    #[test]
    fn threshold_rises_with_electrode_distance() {
        let (t05, t1, t2) = (threshold_ua(0.5), threshold_ua(1.0), threshold_ua(2.0));
        let slope = (t2.ln() - t05.ln()) / (2.0_f64 / 0.5).ln();
        eprintln!("thresholds µA: 0.5mm={t05:.1} 1mm={t1:.1} 2mm={t2:.1} (exponent≈{slope:.2})");
        assert!(
            t05 < t1 && t1 < t2,
            "recruitment threshold must rise with electrode distance: {t05} {t1} {t2}"
        );
    }

    #[test]
    fn electrode_recruits_above_threshold_not_below() {
        let scene = Scene::single_axon(1.0);
        let below = stimulate(&scene, -20.0).expect("solve");
        let above = stimulate(&scene, -400.0).expect("solve");
        assert!(
            !below.any_fired(),
            "weak cathodic stimulus must not recruit"
        );
        assert!(
            above.any_fired(),
            "strong cathodic stimulus must recruit the near axon"
        );
    }

    #[test]
    fn recruitment_curve_is_monotonic() {
        let scene = Scene::bundle(10, 2.0);
        let curve = recruitment_curve(&scene, &[10.0, 50.0, 200.0, 1000.0, 3000.0]).expect("solve");
        let fracs: Vec<f64> = curve.iter().map(|&(_, f)| f).collect();
        for w in fracs.windows(2) {
            assert!(
                w[1] >= w[0],
                "recruitment must not decrease with current: {fracs:?}"
            );
        }
        assert!(
            *fracs.last().unwrap() > 0.0,
            "strong stimulus should recruit someone: {fracs:?}"
        );
    }
}
