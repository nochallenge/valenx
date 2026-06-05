//! The human-readable wind-tunnel run report.
//!
//! An [`AeroResult`] is a big struct; an [`AeroReport`] is the
//! one-page summary a human (or an LLM relaying to a human) actually
//! reads — the headline coefficients, the case conditions, the
//! convergence verdict, and the honest caveats that qualify the
//! numbers.

use crate::api::AeroResult;
use crate::compressible::FlowRegime;
use crate::forces::turbulence_note;

/// A structured, human-readable summary of a wind-tunnel run.
#[derive(Clone, Debug)]
pub struct AeroReport {
    /// The drag coefficient.
    pub cd: f64,
    /// The lift coefficient.
    pub cl: f64,
    /// The side-force coefficient.
    pub cs: f64,
    /// The pitch-moment coefficient.
    pub cm: f64,
    /// The pressure-drag fraction of the total drag (0–1).
    pub pressure_drag_fraction: f64,
    /// The drag area `Cd·A` (m²).
    pub drag_area: f64,
    /// The free-stream dynamic pressure `q∞ = ½·ρ·U∞²` (Pa) of the run.
    pub dynamic_pressure: f64,
    /// The reference area `A` (m²) the coefficients normalise against — the
    /// body's frontal silhouette.
    pub reference_area: f64,
    /// The free-stream Reynolds number.
    pub reynolds_number: f64,
    /// The free-stream Mach number.
    pub mach_number: f64,
    /// Whether the solve converged.
    pub converged: bool,
    /// The number of solver iterations performed.
    pub iterations: usize,
    /// The final residual.
    pub residual: f64,
    /// The mean dimensionless wall distance `y+` over the body.
    pub y_plus_mean: f64,
    /// The peak velocity deficit in the wake.
    pub wake_peak_deficit: f64,
    /// Honest caveat lines qualifying the result.
    pub caveats: Vec<String>,
}

/// The Prandtl–Glauert compressibility factor `1/√(1 − M²)` at Mach `mach` — the
/// subsonic correction that scales a thin-body's incompressible aerodynamic
/// coefficients toward their compressible values (`Cl ≈ Cl_incompressible ·
/// factor`). It is `1` at `M = 0` and grows toward `M = 1`, where the linearised
/// theory diverges; valid only for `0 ≤ M < 1`. Returns `0` outside that range
/// (sonic / supersonic / non-finite), where the correction does not apply.
pub fn prandtl_glauert_factor(mach: f64) -> f64 {
    if mach.is_finite() && (0.0..1.0).contains(&mach) {
        1.0 / (1.0 - mach * mach).sqrt()
    } else {
        0.0
    }
}

/// The finite-wing (3-D) lift-curve slope `a = a₀ / (1 + a₀/(π·e·AR))` (per
/// radian) from Prandtl's lifting-line theory — how a wing of finite aspect
/// ratio `aspect_ratio` and span efficiency `span_efficiency` develops a
/// *gentler* lift slope than its 2-D airfoil section `section_slope_per_rad`
/// (typically ≈ `2π`/rad). The downwash induced by the trailing vortices tilts
/// the local flow and cuts the effective angle of attack, so a given incidence
/// makes less lift. This is the lift-side companion to the induced drag
/// (`crate::sweep::PolarCurve::induced_drag_factor`) — both scale with the same
/// `1/(π·e·AR)` finite-span group. As `AR → ∞` the downwash vanishes and the
/// 2-D slope `a₀` is recovered. Returns `0` for non-physical inputs
/// (`a₀ < 0`, `AR ≤ 0`, `e ≤ 0`, or any non-finite), where the relation does
/// not apply.
pub fn finite_wing_lift_slope(
    section_slope_per_rad: f64,
    aspect_ratio: f64,
    span_efficiency: f64,
) -> f64 {
    if !section_slope_per_rad.is_finite()
        || !aspect_ratio.is_finite()
        || !span_efficiency.is_finite()
        || section_slope_per_rad < 0.0
        || aspect_ratio <= 0.0
        || span_efficiency <= 0.0
    {
        return 0.0;
    }
    section_slope_per_rad
        / (1.0 + section_slope_per_rad / (std::f64::consts::PI * span_efficiency * aspect_ratio))
}

impl AeroReport {
    /// Build a report from a completed [`AeroResult`].
    pub fn from_result(result: &AeroResult) -> AeroReport {
        let coeff = &result.coefficients;
        let total_cd = coeff.cd.abs().max(1e-12);
        let pressure_drag_fraction =
            (coeff.cd_pressure / total_cd).clamp(0.0, 1.0);

        let mut caveats = Vec::new();
        // The standing v1 caveat.
        caveats.push(
            "immersed-boundary Cartesian CFD — a real v1, not ANSYS Fluent / \
             STAR-CCM+ parity: no body-fitted meshing, no DES/LES, coefficient \
             accuracy improves with grid resolution"
                .to_string(),
        );
        // Turbulence model note.
        caveats.push(format!(
            "turbulence: {}",
            turbulence_note(result.flow.turbulence.model)
        ));
        // Convergence note.
        if !result.converged {
            caveats.push(format!(
                "the solve did NOT reach the convergence tolerance \
                 (residual {:.2e} after {} iterations) — treat the \
                 coefficients as provisional",
                result.flow.residual, result.flow.iterations
            ));
        }
        // y+ note — a wall function wants y+ in roughly 30–300.
        let yp = result.surface.y_plus_mean;
        if yp > 0.0 && !(1.0..=1000.0).contains(&yp) {
            caveats.push(format!(
                "mean y+ is {yp:.1} — outside the wall-function-friendly band; \
                 the near-wall resolution does not match the high-Re wall \
                 treatment well"
            ));
        }
        // Compressibility note.
        let regime = FlowRegime::classify(result.mach_number);
        if regime != FlowRegime::Incompressible {
            caveats.push(format!("Mach {:.2}: {}", result.mach_number, regime.caveat()));
        }

        AeroReport {
            cd: coeff.cd,
            cl: coeff.cl,
            cs: coeff.cs,
            cm: coeff.cmy,
            pressure_drag_fraction,
            drag_area: result.drag_area(),
            dynamic_pressure: result.tunnel.dynamic_pressure(),
            reference_area: result.tunnel.reference_area,
            reynolds_number: result.reynolds_number,
            mach_number: result.mach_number,
            converged: result.converged,
            iterations: result.flow.iterations,
            residual: result.flow.residual,
            y_plus_mean: result.surface.y_plus_mean,
            wake_peak_deficit: result.wake.peak_deficit(),
            caveats,
        }
    }

    /// The dimensional drag **force** in newtons — `F_D = (Cd·A)·q∞`, the
    /// drag area scaled by the free-stream dynamic pressure. This is the load
    /// the body actually feels, the number behind the dimensionless `cd`.
    pub fn drag_force(&self) -> f64 {
        self.drag_area * self.dynamic_pressure
    }

    /// The aerodynamic lift **force** in newtons — `F_L = Cl·A·q∞` (negative
    /// when the body makes downforce). The dimensional companion to `cl`.
    pub fn lift_force(&self) -> f64 {
        self.cl * self.reference_area * self.dynamic_pressure
    }

    /// The magnitude of the **total** aerodynamic force on the body (N) — the
    /// vector sum of lift, drag, and side force, `√(L² + D² + S²)`. This is the
    /// resultant load the structure must actually react; it always meets or
    /// exceeds the largest single component because the others add in
    /// quadrature.
    pub fn resultant_force(&self) -> f64 {
        let lift = self.lift_force();
        let drag = self.drag_force();
        let side = self.cs * self.reference_area * self.dynamic_pressure;
        (lift * lift + drag * drag + side * side).sqrt()
    }

    /// The lift-to-drag ratio `L/D = Cl / Cd` at this operating point — the
    /// single headline aerodynamic-efficiency number (in unpowered flight, the
    /// glide ratio: horizontal distance travelled per unit height lost). Returns
    /// `0` when the drag is non-positive.
    pub fn lift_to_drag(&self) -> f64 {
        if self.cd > 1e-9 {
            self.cl / self.cd
        } else {
            0.0
        }
    }

    /// The glide angle `γ = atan2(Cd, Cl)` (radians) — the descent slope in
    /// unpowered flight, where `tan γ = D/L = 1/(L/D)`. Defined for any sign of
    /// lift via `atan2`: a lifting body glides shallowly (small `γ`), a
    /// non-lifting or draggy body descends steeply (`γ → π/2` and beyond).
    pub fn glide_angle_rad(&self) -> f64 {
        self.cd.atan2(self.cl)
    }

    /// The Prandtl–Glauert compressibility factor at this run's Mach number —
    /// the subsonic correction that scales incompressible coefficients toward
    /// their compressible values (`Cl ≈ Cl_incompressible · factor`). See the
    /// free [`prandtl_glauert_factor`]; `1` at low speed, growing toward `M = 1`,
    /// `0` once sonic/supersonic (the linearised correction breaks down).
    pub fn prandtl_glauert_factor(&self) -> f64 {
        prandtl_glauert_factor(self.mach_number)
    }

    /// Render the report as a plain-text block — the form a CLI prints
    /// or an LLM relays.
    pub fn to_text(&self) -> String {
        let mut s = String::new();
        s.push_str("=== valenx-aero wind-tunnel report ===\n");
        s.push_str(&format!(
            "  Reynolds number : {:.3e}\n",
            self.reynolds_number
        ));
        s.push_str(&format!("  Mach number     : {:.3}\n", self.mach_number));
        s.push_str(&format!(
            "  P-G factor      : {:.3}  (1/\u{221A}(1\u{2212}M\u{00B2}))\n",
            self.prandtl_glauert_factor()
        ));
        s.push_str(&format!(
            "  converged       : {} ({} iterations, residual {:.2e})\n",
            self.converged, self.iterations, self.residual
        ));
        s.push_str("  --- coefficients ---\n");
        s.push_str(&format!("  drag      Cd : {:+.4}\n", self.cd));
        s.push_str(&format!("  lift      Cl : {:+.4}\n", self.cl));
        s.push_str(&format!("  side      Cs : {:+.4}\n", self.cs));
        s.push_str(&format!("  pitch Cm     : {:+.4}\n", self.cm));
        s.push_str(&format!("  lift/drag L/D: {:+.3}\n", self.lift_to_drag()));
        s.push_str(&format!(
            "  glide angle  : {:.2} deg\n",
            self.glide_angle_rad().to_degrees()
        ));
        s.push_str(&format!(
            "  drag area CdA: {:.4} m^2\n",
            self.drag_area
        ));
        s.push_str(&format!(
            "  dynamic press: {:.1} Pa  (q_inf)\n",
            self.dynamic_pressure
        ));
        s.push_str(&format!(
            "  drag force   : {:.1} N\n",
            self.drag_force()
        ));
        s.push_str(&format!(
            "  ref area A   : {:.4} m^2\n",
            self.reference_area
        ));
        s.push_str(&format!(
            "  lift force   : {:.1} N\n",
            self.lift_force()
        ));
        s.push_str(&format!(
            "  resultant F  : {:.1} N\n",
            self.resultant_force()
        ));
        s.push_str(&format!(
            "  pressure drag: {:.0}% of total ({}% friction)\n",
            100.0 * self.pressure_drag_fraction,
            (100.0 * (1.0 - self.pressure_drag_fraction)).round() as i64
        ));
        s.push_str("  --- flow ---\n");
        s.push_str(&format!("  mean y+         : {:.1}\n", self.y_plus_mean));
        s.push_str(&format!(
            "  wake peak deficit: {:.0}%\n",
            100.0 * self.wake_peak_deficit
        ));
        s.push_str("  --- caveats ---\n");
        for c in &self.caveats {
            s.push_str(&format!("  * {c}\n"));
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api::{run_windtunnel, AeroRequest};
    use crate::domain::TunnelSizing;
    use crate::geometry::box_body;
    use crate::turbulence::TurbulenceModel;
    use nalgebra::Vector3;

    /// A coarse but real grid for the report-formatting tests — they
    /// assert report *structure* (the headline numbers are copied
    /// through, the caveats are present, the text contains the labels),
    /// none of which needs a fine mesh. A coarse grid keeps the real
    /// end-to-end solve fast.
    fn coarse() -> TunnelSizing {
        TunnelSizing {
            cells_across_body: 4,
            max_cells: 40_000,
            ..TunnelSizing::default()
        }
    }

    #[test]
    fn report_summarises_a_run() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // The headline numbers are carried over verbatim.
        assert_eq!(report.cd, result.coefficients.cd);
        assert_eq!(report.cl, result.coefficients.cl);
        // The pressure-drag fraction is a fraction.
        assert!((0.0..=1.0).contains(&report.pressure_drag_fraction));
        // There is always at least the standing v1 caveat + the
        // turbulence note.
        assert!(report.caveats.len() >= 2);
        assert!(report
            .caveats
            .iter()
            .any(|c| c.contains("not ANSYS Fluent")));
    }

    #[test]
    fn report_text_contains_the_coefficients() {
        let body = box_body(Vector3::zeros(), Vector3::new(1.0, 1.0, 1.0));
        let req = AeroRequest::new(15.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(20);
        let result = run_windtunnel(&body, &req).unwrap();
        let text = AeroReport::from_result(&result).to_text();
        assert!(text.contains("wind-tunnel report"));
        assert!(text.contains("drag      Cd"));
        assert!(text.contains("Reynolds number"));
        assert!(text.contains("caveats"));
    }

    #[test]
    fn non_converged_run_gets_a_caveat() {
        // A 1-iteration run will not converge — the report must say so.
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_max_iterations(1);
        let result = run_windtunnel(&body, &req).unwrap();
        if !result.converged {
            let report = AeroReport::from_result(&result);
            assert!(
                report.caveats.iter().any(|c| c.contains("did NOT reach")),
                "a non-converged run should be flagged"
            );
        }
    }

    #[test]
    fn report_carries_dynamic_pressure_and_drag_force() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // q∞ is copied straight from the tunnel's free-stream conditions.
        assert_eq!(report.dynamic_pressure, result.tunnel.dynamic_pressure());
        assert!(report.dynamic_pressure > 0.0, "moving air has positive q");
        // Drag force is the drag area scaled by q: F_D = (Cd·A)·q∞.
        let expected = report.drag_area * report.dynamic_pressure;
        assert!((report.drag_force() - expected).abs() < 1e-9);
        // It also surfaces in the text dump.
        assert!(report.to_text().contains("drag force"));
    }

    #[test]
    fn report_carries_reference_area_and_lift_force() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // Reference area is the tunnel's, and drag_area = Cd·A stays consistent.
        assert_eq!(report.reference_area, result.tunnel.reference_area);
        assert!(report.reference_area > 0.0, "a real body has frontal area");
        assert!((report.drag_area - report.cd * report.reference_area).abs() < 1e-9);
        // Lift force is Cl·A·q (definitional), and it surfaces in the text dump.
        let expected = report.cl * report.reference_area * report.dynamic_pressure;
        assert!((report.lift_force() - expected).abs() < 1e-9);
        assert!(report.to_text().contains("lift force"));
    }

    #[test]
    fn resultant_force_is_the_quadrature_sum_of_lift_drag_side() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // Resultant = √(L² + D² + S²), the quadrature sum of the three forces.
        let side = report.cs * report.reference_area * report.dynamic_pressure;
        let expected =
            (report.lift_force().powi(2) + report.drag_force().powi(2) + side * side).sqrt();
        assert!((report.resultant_force() - expected).abs() < 1e-9);
        // It meets or exceeds every single component (quadrature, never smaller).
        assert!(report.resultant_force() >= report.drag_force().abs() - 1e-9);
        assert!(report.resultant_force() >= report.lift_force().abs() - 1e-9);
        // And it surfaces in the text dump.
        assert!(report.to_text().contains("resultant F"));
    }

    #[test]
    fn lift_to_drag_is_the_coefficient_ratio() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // L/D · Cd = Cl identically (the definitional ratio), and it is finite.
        assert!((report.lift_to_drag() * report.cd - report.cl).abs() < 1e-9);
        assert!(report.lift_to_drag().is_finite());
        // It surfaces in the text dump.
        assert!(report.to_text().contains("L/D"));
    }

    #[test]
    fn glide_angle_is_the_arctangent_of_drag_over_lift() {
        let body = box_body(Vector3::zeros(), Vector3::new(2.0, 1.0, 1.0));
        let req = AeroRequest::new(20.0)
            .with_turbulence(TurbulenceModel::KEpsilon)
            .with_sizing(coarse())
            .with_max_iterations(30);
        let result = run_windtunnel(&body, &req).unwrap();
        let report = AeroReport::from_result(&result);
        // γ = atan2(Cd, Cl): a body with positive drag descends at an angle
        // strictly between 0 and π.
        let g = report.glide_angle_rad();
        assert!(g > 0.0 && g < std::f64::consts::PI && g.is_finite(), "γ {g}");
        // For a lifting body, tan γ = D/L = 1/(L/D).
        if report.cl > 0.0 && report.lift_to_drag() > 0.0 {
            assert!((g.tan() - 1.0 / report.lift_to_drag()).abs() < 1e-9, "tan γ = 1/(L/D)");
        }
        // It surfaces in the text dump.
        assert!(report.to_text().contains("glide angle"));
    }

    #[test]
    fn prandtl_glauert_factor_matches_textbook_values() {
        // M = 0 → no compressibility correction.
        assert!((prandtl_glauert_factor(0.0) - 1.0).abs() < 1e-12);
        // M = 0.6 → 1/√(1−0.36) = 1/0.8 = 1.25.
        assert!((prandtl_glauert_factor(0.6) - 1.25).abs() < 1e-12);
        // M = 0.8 → 1/√(1−0.64) = 1/0.6 ≈ 1.6667.
        assert!((prandtl_glauert_factor(0.8) - 1.0 / 0.6).abs() < 1e-12);
        // It rises monotonically through the subsonic range.
        assert!(prandtl_glauert_factor(0.7) > prandtl_glauert_factor(0.3));
        // Sonic / supersonic / out-of-range → 0 (the correction breaks down).
        assert_eq!(prandtl_glauert_factor(1.0), 0.0);
        assert_eq!(prandtl_glauert_factor(1.5), 0.0);
        assert_eq!(prandtl_glauert_factor(-0.1), 0.0);
    }

    #[test]
    fn finite_wing_lift_slope_reduces_below_the_section_value() {
        use std::f64::consts::PI;
        let a0 = 2.0 * PI; // thin-airfoil section slope, per radian
        // AR = 6, e = 1: a = 2π / (1 + 2π/(π·6)) = 2π / (1 + 1/3) = 2π·0.75.
        let a = finite_wing_lift_slope(a0, 6.0, 1.0);
        assert!((a - 2.0 * PI * 0.75).abs() < 1e-9, "AR=6 slope {a}");
        // A finite wing is always gentler than its 2-D section.
        assert!(a < a0, "finite-wing slope must be < section slope");

        // As AR → ∞ the downwash vanishes and the 2-D slope is recovered.
        let a_inf = finite_wing_lift_slope(a0, 1.0e6, 1.0);
        assert!((a_inf - a0).abs() < 1e-3, "AR→∞ should recover a0, got {a_inf}");

        // Monotonic: a higher-aspect-ratio wing has a steeper slope.
        assert!(finite_wing_lift_slope(a0, 12.0, 1.0) > finite_wing_lift_slope(a0, 6.0, 1.0));

        // Non-physical inputs → 0 (the relation does not apply).
        assert_eq!(finite_wing_lift_slope(a0, 0.0, 1.0), 0.0);
        assert_eq!(finite_wing_lift_slope(a0, 6.0, 0.0), 0.0);
        assert_eq!(finite_wing_lift_slope(-1.0, 6.0, 1.0), 0.0);
        assert_eq!(finite_wing_lift_slope(f64::NAN, 6.0, 1.0), 0.0);
    }
}
