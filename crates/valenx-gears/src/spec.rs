//! [`GearSpec`] + [`GearKind`].

use serde::{Deserialize, Serialize};

/// Type of gear to generate.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum GearKind {
    /// Parallel-axis straight-tooth (spur).
    Spur,
    /// Parallel-axis with helix angle.
    Helical,
    /// Intersecting-axis truncated-cone bevel.
    Bevel,
    /// Crossed-axis worm (helical thread + worm gear pinion).
    Worm,
}

impl GearKind {
    /// Short UI label.
    pub fn label(self) -> &'static str {
        match self {
            Self::Spur => "Spur",
            Self::Helical => "Helical",
            Self::Bevel => "Bevel",
            Self::Worm => "Worm",
        }
    }
}

/// Parametric description of a single gear.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GearSpec {
    /// Which family.
    pub kind: GearKind,
    /// Module (mm). Pitch diameter = module × teeth.
    pub module_mm: f64,
    /// Tooth count.
    pub teeth: u32,
    /// Pressure angle, degrees. Standard = 20°.
    pub pressure_angle_deg: f64,
    /// Helix angle, degrees. 0 for spur. ~20-30 for helical.
    pub helix_angle_deg: f64,
    /// Face width, mm.
    pub face_width_mm: f64,
}

impl GearSpec {
    /// Convenience: a standard 1-module, 20° spur gear.
    pub fn standard_spur(teeth: u32) -> Self {
        Self {
            kind: GearKind::Spur,
            module_mm: 1.0,
            teeth,
            pressure_angle_deg: 20.0,
            helix_angle_deg: 0.0,
            face_width_mm: 10.0,
        }
    }

    /// Pitch (reference) circle diameter — `module × teeth`.
    pub fn pitch_diameter_mm(&self) -> f64 {
        self.module_mm * self.teeth as f64
    }

    /// Base circle diameter — `pitch × cos(pressure_angle)`. The
    /// involute curve is generated from this circle.
    pub fn base_diameter_mm(&self) -> f64 {
        self.pitch_diameter_mm() * self.pressure_angle_deg.to_radians().cos()
    }

    /// Addendum diameter — pitch + 2 × module (standard).
    pub fn addendum_diameter_mm(&self) -> f64 {
        self.pitch_diameter_mm() + 2.0 * self.module_mm
    }

    /// Dedendum diameter — pitch − 2.5 × module (clearance 0.25 m).
    pub fn dedendum_diameter_mm(&self) -> f64 {
        (self.pitch_diameter_mm() - 2.5 * self.module_mm).max(0.0)
    }
}

impl Default for GearSpec {
    fn default() -> Self {
        Self::standard_spur(20)
    }
}

/// Circular pitch `p = π·m` (mm) — the arc length between corresponding points on adjacent
/// teeth, measured along the pitch circle, for a gear of module `module_mm`. Returns `0.0` for
/// a non-positive or non-finite module.
pub fn circular_pitch_mm(module_mm: f64) -> f64 {
    if !module_mm.is_finite() || module_mm <= 0.0 {
        return 0.0;
    }
    module_mm * std::f64::consts::PI
}

/// Gear ratio `N_driven / N_driver` — the dimensionless speed-reduction (equivalently
/// torque-multiplication) factor of a meshing gear pair. A ratio > 1 reduces speed and
/// multiplies torque; < 1 is an overdrive. Returns `0.0` when `driver_teeth` is 0.
pub fn gear_ratio(driven_teeth: u32, driver_teeth: u32) -> f64 {
    if driver_teeth == 0 {
        return 0.0;
    }
    driven_teeth as f64 / driver_teeth as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circular_pitch_is_pi_times_module() {
        // p = π·m: module 2 → 2π ≈ 6.283185.
        assert!((circular_pitch_mm(2.0) - 2.0 * std::f64::consts::PI).abs() < 1e-12);
        // Linear in module.
        assert!((circular_pitch_mm(4.0) - 2.0 * circular_pitch_mm(2.0)).abs() < 1e-12);
        // Guards: non-positive or non-finite → 0.
        assert_eq!(circular_pitch_mm(0.0), 0.0);
        assert_eq!(circular_pitch_mm(-1.5), 0.0);
        assert_eq!(circular_pitch_mm(f64::NAN), 0.0);
    }

    #[test]
    fn gear_ratio_is_driven_over_driver() {
        // 40:10 → 4:1 reduction; 10:40 → 0.25 overdrive; 20:20 → 1.0.
        assert!((gear_ratio(40, 10) - 4.0).abs() < 1e-12);
        assert!((gear_ratio(10, 40) - 0.25).abs() < 1e-12);
        assert!((gear_ratio(20, 20) - 1.0).abs() < 1e-12);
        // Reciprocal: ratio(a,b)·ratio(b,a) = 1 for nonzero.
        assert!((gear_ratio(40, 10) * gear_ratio(10, 40) - 1.0).abs() < 1e-12);
        // Guard: zero driver teeth → 0.0.
        assert_eq!(gear_ratio(50, 0), 0.0);
    }

    #[test]
    fn base_pitch_equals_pi_m_cos_alpha() {
        // Base pitch p_b = π·m·cos(α) — the spacing of involute flanks
        // measured on the base circle. For module m=2 mm, α=20°:
        //   p_b = π·2·cos(20°) = 5.90430 mm.
        // The crate exposes no base-pitch fn, but base pitch = (base-circle
        // circumference)/teeth = π·base_diameter/N, and base_diameter and the
        // circular pitch p = π·m are both public. Compute p_b two equivalent
        // ways from the real API and pin BOTH to the exact closed form (tol
        // 1e-3 mm, a tight bound for this exact trig form).
        let mut g = GearSpec::standard_spur(20);
        g.module_mm = 2.0; // 20° pressure angle is the standard_spur default
        let expected = std::f64::consts::PI * 2.0 * 20.0_f64.to_radians().cos(); // 5.90430 mm
        assert!(
            (expected - 5.904_30).abs() < 1e-4,
            "ground-truth check: {expected}"
        );

        // Route 1: via base circle — p_b = π·D_base / N.
        let pb_via_base = std::f64::consts::PI * g.base_diameter_mm() / g.teeth as f64;
        assert!(
            (pb_via_base - expected).abs() < 1e-3,
            "base pitch via base circle {pb_via_base} mm vs closed form {expected} mm"
        );

        // Route 2: via circular pitch — p_b = p·cos(α) = (π·m)·cos(α).
        let pb_via_cp = circular_pitch_mm(g.module_mm) * g.pressure_angle_deg.to_radians().cos();
        assert!(
            (pb_via_cp - expected).abs() < 1e-3,
            "base pitch via circular pitch {pb_via_cp} mm vs closed form {expected} mm"
        );

        // Base pitch is strictly less than circular pitch (cos α < 1).
        assert!(pb_via_cp < circular_pitch_mm(g.module_mm));
    }

    #[test]
    fn diameters_follow_standard_gear_relations() {
        // module = 2 (NOT the default 1): with module 1, `module·teeth`,
        // `2·module`, and `2.5·module` collapse to teeth, 2, 2.5 — so a
        // regression dropping the module factor would pass any module-1 test.
        let mut g = GearSpec::standard_spur(20);
        g.module_mm = 2.0;
        let pitch = g.pitch_diameter_mm();
        assert!((pitch - 40.0).abs() < 1e-12); // module·teeth = 2·20

        // Base circle = pitch·cos(pressure_angle); 20° here.
        let base = g.base_diameter_mm();
        assert!((base - 40.0 * 20.0_f64.to_radians().cos()).abs() < 1e-12);

        // Addendum diameter − pitch = 2·module; pitch − dedendum = 2.5·module.
        assert!((g.addendum_diameter_mm() - pitch - 2.0 * g.module_mm).abs() < 1e-12);
        assert!((pitch - g.dedendum_diameter_mm() - 2.5 * g.module_mm).abs() < 1e-12);

        // Standard diameter ordering: dedendum < base < pitch < addendum.
        assert!(g.dedendum_diameter_mm() < base);
        assert!(base < pitch);
        assert!(pitch < g.addendum_diameter_mm());

        // Pitch is linear in BOTH teeth and module (reconstruct, not clone).
        let mut g_teeth = GearSpec::standard_spur(40);
        g_teeth.module_mm = 2.0; // teeth doubled → pitch doubles
        assert!((g_teeth.pitch_diameter_mm() - 2.0 * pitch).abs() < 1e-12);
        let mut g_mod = GearSpec::standard_spur(20);
        g_mod.module_mm = 4.0; // module doubled → pitch doubles
        assert!((g_mod.pitch_diameter_mm() - 2.0 * pitch).abs() < 1e-12);

        // Dedendum clamps to 0 for a tiny gear where pitch < 2.5·module
        // (1 tooth, module 1 → pitch 1 < 2.5 → the `.max(0.0)` engages).
        assert_eq!(GearSpec::standard_spur(1).dedendum_diameter_mm(), 0.0);
    }
}
