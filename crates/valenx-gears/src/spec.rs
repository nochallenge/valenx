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

    /// Base pitch `p_b = π·m·cos(α)` (mm) — the spacing of involute tooth
    /// flanks measured along the base circle, equal to the circular pitch
    /// `π·m` reduced by `cos(pressure angle)`. It is the pitch that governs
    /// how successive teeth hand off contact, and so sets the
    /// [`contact_ratio`]. Returns `0.0` for a non-positive / non-finite
    /// module (inherited from [`circular_pitch_mm`]).
    pub fn base_pitch_mm(&self) -> f64 {
        circular_pitch_mm(self.module_mm) * self.pressure_angle_deg.to_radians().cos()
    }

    /// Lewis tooth-root **bending stress** `σ` (MPa) for this gear,
    /// carrying a tangential transmitted load `w_t_n` (newtons) on a tooth
    /// whose Lewis form factor is `lewis_form_factor_y`.
    ///
    /// Convenience wrapper over [`lewis_bending_stress_mpa`] that uses this
    /// gear's own module and face width:
    ///
    /// ```text
    /// σ = W_t / (F · m · Y)
    /// ```
    ///
    /// See [`lewis_bending_stress_mpa`] for the full definition, units and
    /// the Shigley reference. Returns `0.0` for non-positive / non-finite
    /// inputs.
    pub fn lewis_bending_stress_mpa(&self, w_t_n: f64, lewis_form_factor_y: f64) -> f64 {
        lewis_bending_stress_mpa(
            w_t_n,
            self.face_width_mm,
            self.module_mm,
            lewis_form_factor_y,
        )
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

/// **Lewis bending stress** `σ` (MPa) at the root of a gear tooth,
/// modelled as a cantilever beam of uniform strength loaded at its tip
/// by the tangential transmitted force.
///
/// ```text
/// σ = W_t / (F · m · Y)
/// ```
///
/// with, in a consistent SI-millimetre system,
///
/// * `w_t_n` — tangential transmitted load `W_t` (N),
/// * `face_width_mm` — face width `F` (mm),
/// * `module_mm` — module `m` (mm),
/// * `lewis_form_factor_y` — the dimensionless Lewis form factor `Y`,
///   read from a table by tooth count and pressure angle (e.g. Shigley's
///   *Mechanical Engineering Design*, Table 14-2; for 20 teeth at a 20°
///   full-depth tooth, `Y ≈ 0.322`).
///
/// Because `N·mm⁻² = MPa`, the result is in megapascals directly. This is
/// the original Lewis equation (Shigley, *Mechanical Engineering Design*,
/// "The Lewis Bending Equation") — the bare geometric/strength estimate
/// **without** the AGMA velocity, size, load-distribution or geometry
/// (`J`) refinement factors, so it is a first-order screening value, not a
/// rated allowable. Returns `0.0` for any non-positive or non-finite
/// argument (the crate's invalid-input sentinel).
pub fn lewis_bending_stress_mpa(
    w_t_n: f64,
    face_width_mm: f64,
    module_mm: f64,
    lewis_form_factor_y: f64,
) -> f64 {
    if ![w_t_n, face_width_mm, module_mm, lewis_form_factor_y]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
    {
        return 0.0;
    }
    w_t_n / (face_width_mm * module_mm * lewis_form_factor_y)
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

/// Transverse **contact ratio** `mₚ` of a meshing spur-gear pair at standard
/// center distance — the average number of tooth pairs sharing the load.
///
/// ```text
/// mₚ = [ √(ra1² − rb1²) + √(ra2² − rb2²) − C·sin(α) ] / p_b
/// ```
///
/// where `ra` is the addendum (outside) radius, `rb` the base radius, `C` the
/// center distance `(d1 + d2)/2`, `α` the pressure angle and `p_b` the
/// [base pitch](GearSpec::base_pitch_mm). The numerator is the length of the
/// active line of action; dividing by the base pitch counts how many tooth
/// engagements span it. A pair **must** have `mₚ > 1` for continuous
/// transmission (always at least one tooth in contact); typical standard
/// spur pairs sit around `1.4–1.8`.
///
/// Two involute gears mesh only at a **common module and pressure angle**;
/// this returns `0.0` (the crate's invalid-input sentinel, as with
/// [`circular_pitch_mm`] / [`gear_ratio`]) when they differ or when the
/// module is non-positive / non-finite.
pub fn contact_ratio(pinion: &GearSpec, gear: &GearSpec) -> f64 {
    let m = pinion.module_mm;
    if !m.is_finite()
        || m <= 0.0
        || (pinion.module_mm - gear.module_mm).abs() > 1e-9
        || (pinion.pressure_angle_deg - gear.pressure_angle_deg).abs() > 1e-9
    {
        return 0.0;
    }
    let pb = pinion.base_pitch_mm();
    if pb <= 0.0 {
        return 0.0;
    }
    let phi = pinion.pressure_angle_deg.to_radians();
    let (ra1, rb1) = (
        pinion.addendum_diameter_mm() / 2.0,
        pinion.base_diameter_mm() / 2.0,
    );
    let (ra2, rb2) = (
        gear.addendum_diameter_mm() / 2.0,
        gear.base_diameter_mm() / 2.0,
    );
    // Standard center distance C = (d1 + d2)/2.
    let c = (pinion.pitch_diameter_mm() + gear.pitch_diameter_mm()) / 2.0;
    let active_line = (ra1 * ra1 - rb1 * rb1).max(0.0).sqrt()
        + (ra2 * ra2 - rb2 * rb2).max(0.0).sqrt()
        - c * phi.sin();
    active_line / pb
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
    fn lewis_bending_stress_matches_worked_example() {
        // Worked example (Lewis equation σ = W_t/(F·m·Y), SI form):
        //   module m = 2 mm, face width F = 20 mm, tangential load
        //   W_t = 500 N, Lewis form factor Y = 0.35 (20-tooth gear).
        //   σ = 500 / (20 · 2 · 0.35) = 35.714 MPa.
        // Reference: Lewis bending equation, Shigley's *Mechanical
        // Engineering Design*; this exact worked example (35.7 MPa) is
        // reproduced by the FIRGELLI Lewis gear-tooth-strength calculator
        // <https://www.firgelliauto.com/blogs/engineering-calculators/gear-tooth-strength-calculator-lewis-formula>.
        let sigma = lewis_bending_stress_mpa(500.0, 20.0, 2.0, 0.35);
        assert!(
            (sigma - 35.714_285_7).abs() < 1e-4,
            "Lewis stress {sigma} MPa, expected 35.714 MPa"
        );

        // Same load/geometry but with the Shigley Table 14-2 Lewis form
        // factor for a 20-tooth, 20° full-depth tooth, Y = 0.322:
        //   σ = 500 / (20 · 2 · 0.322) = 38.820 MPa.
        // (Shigley's *Mechanical Engineering Design*, Table 14-2.)
        let sigma_shigley = lewis_bending_stress_mpa(500.0, 20.0, 2.0, 0.322);
        assert!(
            (sigma_shigley - 38.819_875).abs() < 1e-4,
            "Lewis stress {sigma_shigley} MPa, expected 38.820 MPa"
        );

        // The GearSpec convenience method uses the gear's own module and
        // face width and must agree with the free function.
        let mut g = GearSpec::standard_spur(20);
        g.module_mm = 2.0;
        g.face_width_mm = 20.0;
        assert!((g.lewis_bending_stress_mpa(500.0, 0.322) - sigma_shigley).abs() < 1e-12);

        // Physical scalings: stress is inversely proportional to face
        // width, module and Y, and linear in load.
        let base = lewis_bending_stress_mpa(500.0, 20.0, 2.0, 0.322);
        assert!((lewis_bending_stress_mpa(1000.0, 20.0, 2.0, 0.322) - 2.0 * base).abs() < 1e-9);
        assert!((lewis_bending_stress_mpa(500.0, 40.0, 2.0, 0.322) - 0.5 * base).abs() < 1e-9);

        // Guards: any non-positive / non-finite argument → 0.0 sentinel.
        assert_eq!(lewis_bending_stress_mpa(0.0, 20.0, 2.0, 0.322), 0.0);
        assert_eq!(lewis_bending_stress_mpa(500.0, 0.0, 2.0, 0.322), 0.0);
        assert_eq!(lewis_bending_stress_mpa(500.0, 20.0, -2.0, 0.322), 0.0);
        assert_eq!(lewis_bending_stress_mpa(500.0, 20.0, 2.0, f64::NAN), 0.0);
    }

    #[test]
    fn base_pitch_equals_pi_m_cos_alpha() {
        // Base pitch p_b = π·m·cos(α) — the spacing of involute flanks
        // measured on the base circle. For module m=2 mm, α=20°:
        //   p_b = π·2·cos(20°) = 5.90430 mm.
        // Base pitch is exposed by `base_pitch_mm()`, and can also be derived
        // two other ways: base pitch = (base-circle circumference)/teeth =
        // π·base_diameter/N, and base pitch = circular pitch × cos(α). Compute
        // p_b all three ways from the real API and pin each to the exact
        // closed form (tol 1e-3 mm, a tight bound for this exact trig form).
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

        // Route 3: the dedicated base_pitch_mm() accessor agrees with both.
        assert!(
            (g.base_pitch_mm() - expected).abs() < 1e-3,
            "base_pitch_mm() {} mm vs closed form {expected} mm",
            g.base_pitch_mm()
        );

        // Base pitch is strictly less than circular pitch (cos α < 1).
        assert!(pb_via_cp < circular_pitch_mm(g.module_mm));
    }

    #[test]
    fn contact_ratio_of_two_standard_20_tooth_gears() {
        // Two identical 20-tooth, module-1, 20° spur gears at standard center
        // distance. Hand computation:
        //   r=10, rb=10·cos20°=9.39693, ra=11, C=20, p_b=π·cos20°=2.952155
        //   mp = [2·√(11²−9.39693²) − 20·sin20°] / 2.952155
        //      = [2·5.71819 − 6.84040] / 2.952155 = 1.55681
        let g = GearSpec::standard_spur(20); // module 1, 20° by default
        let mp = contact_ratio(&g, &g);
        assert!((mp - 1.556_81).abs() < 1e-4, "contact ratio {mp}");
        // A meshing pair must keep at least one tooth pair engaged.
        assert!(mp > 1.0);
    }

    #[test]
    fn contact_ratio_of_dissimilar_pair_20_by_40() {
        // A real reduction pair: 20-tooth pinion meshing a 40-tooth gear,
        // module 1, 20° full-depth, standard center distance. Hand
        // computation (radii in mm):
        //   pinion: r=10, rb=10·cos20°=9.39693, ra=11
        //   gear:   r=20, rb=20·cos20°=18.79385, ra=21
        //   C=30, p_b=π·cos20°=2.95213
        //   active line = √(11²−9.39693²) + √(21²−18.79385²) − 30·sin20°
        //               = 5.71819 + 9.36950 − 10.26060 = 4.82728
        //   m_c = 4.82728 / 2.95213 = 1.63519
        let pinion = GearSpec::standard_spur(20);
        let gear = GearSpec::standard_spur(40);
        let mp = contact_ratio(&pinion, &gear);
        assert!((mp - 1.635_19).abs() < 1e-4, "contact ratio {mp}");
        // Order-independent for a meshing pair.
        assert!((contact_ratio(&gear, &pinion) - mp).abs() < 1e-12);
        // Sits between the identical-20T (≈1.557) and identical-40T values,
        // and stays above 1 for continuous transmission.
        assert!(mp > contact_ratio(&pinion, &pinion));
        assert!(mp > 1.0);
    }

    #[test]
    fn contact_ratio_rises_with_tooth_count_and_rejects_mismatch() {
        // More teeth -> longer line of action -> higher contact ratio.
        let small = GearSpec::standard_spur(18);
        let big = GearSpec::standard_spur(60);
        let mp_small = contact_ratio(&small, &small);
        let mp_big = contact_ratio(&big, &big);
        assert!(mp_big > mp_small, "{mp_big} should exceed {mp_small}");
        assert!(mp_big > 1.0 && mp_small > 1.0);

        // Gears that cannot mesh (different module, or different pressure
        // angle) return the 0.0 sentinel.
        let mut other_module = GearSpec::standard_spur(20);
        other_module.module_mm = 2.0;
        assert_eq!(
            contact_ratio(&GearSpec::standard_spur(20), &other_module),
            0.0
        );
        let mut other_angle = GearSpec::standard_spur(20);
        other_angle.pressure_angle_deg = 25.0;
        assert_eq!(
            contact_ratio(&GearSpec::standard_spur(20), &other_angle),
            0.0
        );
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
