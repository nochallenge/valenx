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

/// The level-flight **stall speed** `V_stall = √(2W / (ρ·S·C_Lmax))` (m/s) — the
/// slowest speed at which a wing of area `wing_area_m2` (m²) can still carry the
/// weight `weight_n` (N) in level flight at air density `air_density` (kg/m³),
/// flying at its maximum lift coefficient `cl_max`. It is exactly the speed at
/// which the available lift `L = ½ρV²S·C_Lmax` equals the weight; any slower and
/// the wing cannot make enough lift and the aircraft stalls. This is the
/// reference speed behind the approach and landing speeds (typically flown at
/// `1.2–1.3·V_stall`). Returns `0` for any non-physical input (non-finite or
/// non-positive).
pub fn stall_speed(weight_n: f64, wing_area_m2: f64, air_density: f64, cl_max: f64) -> f64 {
    if !weight_n.is_finite()
        || weight_n <= 0.0
        || !wing_area_m2.is_finite()
        || wing_area_m2 <= 0.0
        || !air_density.is_finite()
        || air_density <= 0.0
        || !cl_max.is_finite()
        || cl_max <= 0.0
    {
        return 0.0;
    }
    (2.0 * weight_n / (air_density * wing_area_m2 * cl_max)).sqrt()
}

/// The **Breguet range** `R = (V/c)·(L/D)·ln(W₀/W₁)` (m) of a jet in cruise —
/// the classic flight-performance result for how far an aircraft flies burning
/// its fuel. `velocity` `V` (m/s) is the cruise speed, `sfc` `c` (s⁻¹) the
/// thrust-specific fuel consumption (fuel weight per unit thrust per second),
/// `lift_to_drag` `L/D` the aerodynamic efficiency, and `weight_ratio` `W₀/W₁`
/// the start-to-end (takeoff-to-landing) weight ratio — the fuel fraction in
/// log form. It is maximised by cruising at the best-`L/D` point
/// ([`crate::sweep::PolarCurve::best_lift_to_drag_point`]) at the highest
/// `V`/`c`. Returns `0` when no fuel is burned (`W₀/W₁ = 1`) and for any
/// non-physical input (`V`, `c`, or `L/D` non-positive, `W₀/W₁ < 1`, non-finite).
pub fn breguet_range(velocity: f64, sfc: f64, lift_to_drag: f64, weight_ratio: f64) -> f64 {
    if !velocity.is_finite()
        || velocity <= 0.0
        || !sfc.is_finite()
        || sfc <= 0.0
        || !lift_to_drag.is_finite()
        || lift_to_drag <= 0.0
        || !weight_ratio.is_finite()
        || weight_ratio < 1.0
    {
        return 0.0;
    }
    velocity / sfc * lift_to_drag * weight_ratio.ln()
}

/// The **Mach angle** `μ = arcsin(1/M)` (radians) — the half-angle of the Mach
/// cone trailing a body moving at Mach number `mach`, the line along which the
/// pressure disturbances pile into a shock. It is the *supersonic* companion to
/// the subsonic [`prandtl_glauert_factor`]: at `M = 1` the cone is a flat disc
/// (`μ = π/2`), and as the speed rises the cone narrows toward `μ → 0`, sweeping
/// ever more sharply back. There is no Mach cone below the speed of sound, so it
/// returns `0` for subsonic `M < 1` (or non-finite input), exactly as the
/// Prandtl–Glauert factor returns `0` outside its valid range.
pub fn mach_angle(mach: f64) -> f64 {
    if mach.is_finite() && mach >= 1.0 {
        (1.0 / mach).asin()
    } else {
        0.0
    }
}

/// The **induced-drag coefficient** `C_Di = C_L² / (π·e·AR)` of a finite wing
/// (Prandtl lifting-line theory) — the unavoidable "drag-due-to-lift" that comes
/// with making lift at all. A wing of finite aspect ratio `aspect_ratio` `AR`
/// and span efficiency `span_efficiency` `e` sheds trailing vortices whose
/// downwash tilts the lift vector slightly aft, and that backward tilt *is* the
/// induced drag. It is the drag-side companion to [`finite_wing_lift_slope`]:
/// both carry the same `1/(π·e·AR)` finite-span group, here multiplied by the
/// square of the lift coefficient `lift_coefficient` `C_L` (so it quadruples
/// when `C_L` doubles, and is sign-blind — negative-lift downforce induces drag
/// just the same). Elliptical loading (`e = 1`) gives the theoretical minimum
/// `C_Di = C_L²/(π·AR)`, and the term vanishes for an infinite-span wing
/// (`AR → ∞`). Returns `0` at zero lift and for any non-physical input
/// (`AR ≤ 0`, `e ≤ 0`, or non-finite), where the relation does not apply.
pub fn induced_drag_coefficient(
    lift_coefficient: f64,
    aspect_ratio: f64,
    span_efficiency: f64,
) -> f64 {
    if !lift_coefficient.is_finite()
        || !aspect_ratio.is_finite()
        || !span_efficiency.is_finite()
        || aspect_ratio <= 0.0
        || span_efficiency <= 0.0
    {
        return 0.0;
    }
    lift_coefficient * lift_coefficient
        / (std::f64::consts::PI * span_efficiency * aspect_ratio)
}

/// The **maximum lift-to-drag ratio** `(L/D)_max = ½·√(π·e·AR / C_D0)` of a wing
/// with parabolic drag polar `C_D = C_D0 + C_L²/(π·e·AR)` — the single headline
/// aircraft-performance number, the best-glide and best-range/endurance
/// operating point. `cd0` is the zero-lift (parasite) drag coefficient `C_D0`,
/// `aspect_ratio` the wing `AR`, and `span_efficiency` the Oswald factor `e`.
///
/// It is reached at the lift coefficient `C_L* = √(π·e·AR·C_D0)`, the point where
/// the lift-induced drag (see [`induced_drag_coefficient`]) exactly equals the
/// parasite drag `C_D0`; the total drag is then `2·C_D0`, so
/// `(L/D)_max = C_L*/(2·C_D0)`. It improves with span efficiency and aspect ratio
/// (a long, clean wing glides far, `∝ √AR`) and degrades as the parasite drag
/// grows (`∝ 1/√C_D0`). This is the closed-form value from the polar
/// *parameters* — the analytic counterpart to
/// [`crate::sweep::PolarCurve::max_lift_to_drag`], which instead reads the peak
/// off a set of *measured* polar points. Returns `0` for any non-physical input
/// (`C_D0 ≤ 0`, `AR ≤ 0`, `e ≤ 0`, or non-finite).
pub fn max_lift_to_drag_ratio(cd0: f64, aspect_ratio: f64, span_efficiency: f64) -> f64 {
    if !cd0.is_finite()
        || cd0 <= 0.0
        || !aspect_ratio.is_finite()
        || aspect_ratio <= 0.0
        || !span_efficiency.is_finite()
        || span_efficiency <= 0.0
    {
        return 0.0;
    }
    0.5 * (std::f64::consts::PI * span_efficiency * aspect_ratio / cd0).sqrt()
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

    #[test]
    fn stall_speed_is_where_lift_equals_weight_at_max_lift() {
        // A light aircraft: 11 kN weight, 16 m² wing, sea-level air, C_Lmax = 1.5.
        let (w, s, rho, cl_max) = (11_000.0, 16.0, 1.225, 1.5);
        let v = stall_speed(w, s, rho, cl_max);
        assert!((v - 27.36).abs() < 0.1, "stall speed {v} m/s");
        // By construction the lift at V_stall and C_Lmax exactly balances the weight.
        let lift = 0.5 * rho * v * v * s * cl_max;
        assert!((lift - w).abs() < 1e-6, "L={lift} must equal W={w} at the stall");
        // V_stall ∝ 1/√C_Lmax: a higher max lift coefficient lowers the stall speed.
        let v_flapped = stall_speed(w, s, rho, 2.0 * cl_max);
        assert!((v_flapped - v / 2.0_f64.sqrt()).abs() < 1e-9, "∝ 1/√C_Lmax");
        // V_stall ∝ √W: four times the weight doubles the stall speed.
        let v_heavy = stall_speed(4.0 * w, s, rho, cl_max);
        assert!((v_heavy - 2.0 * v).abs() < 1e-9, "∝ √W");
        // Non-physical inputs → 0.
        assert_eq!(stall_speed(0.0, s, rho, cl_max), 0.0);
        assert_eq!(stall_speed(w, s, rho, 0.0), 0.0);
        assert_eq!(stall_speed(w, -1.0, rho, cl_max), 0.0);
        assert_eq!(stall_speed(f64::NAN, s, rho, cl_max), 0.0);
    }

    #[test]
    fn breguet_range_scales_with_speed_efficiency_and_fuel_fraction() {
        // R = (V/c)·(L/D)·ln(W₀/W₁). A long-haul jet: V = 250 m/s, c = 0.6/hr,
        // L/D = 17, W₀/W₁ = 1.5 → ≈ 1.03e7 m (~10,000 km).
        let v = 250.0_f64;
        let c = 0.6 / 3600.0; // 0.6 per hour → per second
        let ld = 17.0_f64;
        let ratio = 1.5_f64;
        let r = breguet_range(v, c, ld, ratio);
        let expected = v / c * ld * ratio.ln();
        assert!((r - expected).abs() < 1e-6, "range {r} vs {expected}");
        // A sensible long-haul figure (within a few thousand km of 10,000).
        assert!((6.0e6..1.4e7).contains(&r), "airliner range {:.0} km", r / 1e3);

        // Burning no fuel (W₀ = W₁) gives zero range.
        assert_eq!(breguet_range(v, c, ld, 1.0), 0.0);
        // R ∝ L/D, ∝ V, and ∝ ln(W₀/W₁) (ratio → ratio² doubles the log).
        assert!((breguet_range(v, c, 2.0 * ld, ratio) - 2.0 * r).abs() < 1e-6, "∝ L/D");
        assert!((breguet_range(2.0 * v, c, ld, ratio) - 2.0 * r).abs() < 1e-6, "∝ V");
        assert!((breguet_range(v, c, ld, ratio * ratio) - 2.0 * r).abs() < 1e-6, "∝ ln(W₀/W₁)");

        // Non-physical inputs → 0 (including an unphysical W₁ > W₀).
        assert_eq!(breguet_range(v, 0.0, ld, ratio), 0.0);
        assert_eq!(breguet_range(-1.0, c, ld, ratio), 0.0);
        assert_eq!(breguet_range(v, c, ld, 0.5), 0.0);
        assert_eq!(breguet_range(v, c, ld, f64::NAN), 0.0);
    }

    #[test]
    fn mach_angle_is_the_mach_cone_half_angle() {
        use std::f64::consts::PI;
        // At M = 1 the Mach cone is a flat disc: μ = 90°.
        assert!((mach_angle(1.0) - PI / 2.0).abs() < 1e-12);
        // M = 2 → arcsin(0.5) = 30°; M = √2 → arcsin(1/√2) = 45°.
        assert!((mach_angle(2.0) - PI / 6.0).abs() < 1e-12, "M=2 → 30°");
        assert!((mach_angle(2.0_f64.sqrt()) - PI / 4.0).abs() < 1e-12, "M=√2 → 45°");
        // The cone narrows as Mach rises, and → 0 at hypersonic speed.
        assert!(mach_angle(5.0) < mach_angle(2.0), "narrows with Mach");
        assert!(mach_angle(1.0e6) < 1e-3, "→ 0 as M → ∞");
        // Subsonic (no Mach cone) and non-finite → 0.
        assert_eq!(mach_angle(0.8), 0.0);
        assert_eq!(mach_angle(0.0), 0.0);
        assert_eq!(mach_angle(f64::NAN), 0.0);
        assert_eq!(mach_angle(f64::INFINITY), 0.0);
    }

    #[test]
    fn induced_drag_coefficient_is_the_lifting_line_drag_due_to_lift() {
        use std::f64::consts::PI;
        // No lift, no induced drag.
        assert_eq!(induced_drag_coefficient(0.0, 8.0, 0.8), 0.0);
        // Worked point: C_L=0.5, AR=8, e=0.8 → 0.25/(π·6.4) ≈ 0.0124339.
        let cdi = induced_drag_coefficient(0.5, 8.0, 0.8);
        assert!((cdi - 0.25 / (PI * 6.4)).abs() < 1e-12, "C_Di {cdi}");
        assert!((cdi - 0.012_433_9).abs() < 1e-6, "≈0.0124339, got {cdi}");
        // Elliptical loading (e=1) is the theoretical minimum C_Di = C_L²/(π·AR),
        // and it is strictly less than the e=0.8 case.
        let cdi_ell = induced_drag_coefficient(0.5, 8.0, 1.0);
        assert!((cdi_ell - 0.25 / (PI * 8.0)).abs() < 1e-12, "elliptical {cdi_ell}");
        assert!(cdi_ell < cdi, "elliptical loading has the least induced drag");
        // Quadratic in C_L: doubling the lift coefficient quadruples induced drag.
        let base = induced_drag_coefficient(0.4, 8.0, 0.8);
        assert!(
            (induced_drag_coefficient(0.8, 8.0, 0.8) - 4.0 * base).abs() < 1e-12,
            "∝ C_L²"
        );
        // Inverse in AR: doubling the aspect ratio halves the induced drag.
        assert!(
            (induced_drag_coefficient(0.5, 16.0, 0.8) - cdi / 2.0).abs() < 1e-12,
            "∝ 1/AR"
        );
        // Sign-blind: negative-lift downforce induces drag just the same (∝ C_L²).
        assert!(
            (induced_drag_coefficient(-0.5, 8.0, 0.8) - cdi).abs() < 1e-12,
            "C_L² is sign-blind"
        );
        // As AR → ∞ the induced drag vanishes (the downwash disappears).
        assert!(induced_drag_coefficient(0.5, 1.0e6, 1.0) < 1e-6, "AR→∞ → 0");
        // Non-physical inputs → 0 (the relation does not apply).
        assert_eq!(induced_drag_coefficient(0.5, 0.0, 0.8), 0.0);
        assert_eq!(induced_drag_coefficient(0.5, 8.0, 0.0), 0.0);
        assert_eq!(induced_drag_coefficient(f64::NAN, 8.0, 0.8), 0.0);
        assert_eq!(induced_drag_coefficient(0.5, f64::INFINITY, 0.8), 0.0);
    }

    #[test]
    fn max_lift_to_drag_ratio_is_the_best_glide_optimum() {
        use std::f64::consts::PI;
        let (cd0, ar, e) = (0.02, 8.0, 0.8);
        let ld_max = max_lift_to_drag_ratio(cd0, ar, e);
        // Closed form ½·√(π·e·AR/C_D0); worked point ≈ 15.85.
        assert!(
            (ld_max - 0.5 * (PI * e * ar / cd0).sqrt()).abs() < 1e-12,
            "closed form"
        );
        assert!((ld_max - 15.85).abs() < 0.05, "≈15.85, got {ld_max}");
        // The optimum is where induced drag equals parasite drag: at
        // C_L* = √(π·e·AR·C_D0), induced_drag_coefficient(C_L*) = C_D0, the total
        // drag is 2·C_D0 and (L/D)_max = C_L*/(2·C_D0).
        let cl_star = (PI * e * ar * cd0).sqrt();
        let cdi = induced_drag_coefficient(cl_star, ar, e);
        assert!((cdi - cd0).abs() < 1e-12, "induced = parasite at the optimum: {cdi}");
        assert!(
            (ld_max - cl_star / (2.0 * cd0)).abs() < 1e-9,
            "(L/D)_max = C_L*/(2·C_D0)"
        );
        // Scaling: ∝ √AR (4× the aspect ratio doubles it) and ∝ 1/√C_D0
        // (4× the parasite drag halves it).
        assert!(
            (max_lift_to_drag_ratio(cd0, 4.0 * ar, e) - 2.0 * ld_max).abs() < 1e-9,
            "∝ √AR"
        );
        assert!(
            (max_lift_to_drag_ratio(4.0 * cd0, ar, e) - ld_max / 2.0).abs() < 1e-9,
            "∝ 1/√C_D0"
        );
        // Non-physical inputs → 0.
        assert_eq!(max_lift_to_drag_ratio(0.0, ar, e), 0.0);
        assert_eq!(max_lift_to_drag_ratio(cd0, 0.0, e), 0.0);
        assert_eq!(max_lift_to_drag_ratio(cd0, ar, 0.0), 0.0);
        assert_eq!(max_lift_to_drag_ratio(f64::NAN, ar, e), 0.0);
    }
}
