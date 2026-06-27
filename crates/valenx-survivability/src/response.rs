//! **Protective structural response** — how a protective element (a wall,
//! plate, or panel reduced to one effective degree of freedom) *survives* a
//! blast pulse.
//!
//! Two complementary tools, both standard in blast-resistant design:
//!
//! 1. [`sdof_response`] — march the equivalent SDOF system under the
//!    [`Friedlander`](valenx_fem::FriedlanderPulse) pressure pulse and read the
//!    **peak deflection** and the **ductility ratio** `μ = x_max / x_yield`
//!    (the demand-to-capacity measure that classifies protective damage). This
//!    reuses the validated Newmark-β transient integrator in
//!    [`valenx_fem`] — we do not reimplement time integration.
//!
//! 2. [`PressureImpulseDiagram`] — the **pressure–impulse (P–I) diagram**, the
//!    canonical protective-design chart. For a chosen damage threshold (an
//!    elastic-energy / deflection limit) it draws the iso-damage curve in
//!    pressure–impulse space: any `(P, i)` *below/left* of the curve is safe,
//!    any point *above/right* exceeds the limit. The curve has two physical
//!    asymptotes — a vertical **impulsive** asymptote (very short loads: only
//!    impulse matters) and a horizontal **quasi-static** asymptote (very long
//!    loads: only peak pressure matters). Pinning those two asymptotes is the
//!    correctness check for the whole diagram.
//!
//! All framing is defensive: these answer "does this protective element survive
//! the design load, and what is the minimum protection that does?" — never how
//! to defeat a structure.

use crate::error::SurvivabilityError;
use serde::{Deserialize, Serialize};
use valenx_fem::{solve_sdof_blast, FriedlanderPulse, TransientControls};

/// Validate a finite, strictly-positive quantity, naming it on error.
fn require_pos(name: &str, v: f64) -> Result<f64, SurvivabilityError> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(SurvivabilityError::InvalidParameter(format!(
            "{name} must be finite and > 0, got {v}"
        )))
    }
}

/// The result of an SDOF protective-element response solve.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SdofResponse {
    /// Peak absolute deflection `x_max` (m) reached over the response history.
    pub peak_deflection_m: f64,
    /// Peak absolute acceleration (m/s²) — feeds the occupant-injury screen.
    pub peak_acceleration_m_s2: f64,
    /// Ductility ratio `μ = x_max / x_yield`. `μ ≤ 1` means the element stayed
    /// elastic (undamaged); `μ > 1` quantifies inelastic demand. `None` if no
    /// yield deflection was supplied.
    pub ductility_ratio: Option<f64>,
    /// Equivalent static deflection `x_st = P_peak / k` (m), the deflection the
    /// peak load would cause if applied slowly — the DAF reference.
    pub static_deflection_m: f64,
    /// Dynamic amplification factor `x_max / x_st`.
    pub dynamic_amplification: f64,
}

/// March an equivalent SDOF protective element (effective mass `m` (kg),
/// damping `c` (N·s/m), stiffness `k` (N/m)) under the Friedlander `pulse`,
/// returning the peak response and ductility.
///
/// `yield_deflection_m`, if `Some`, is the deflection at which the element first
/// yields; it sets the ductility ratio. The `pulse.peak_overpressure` is taken
/// as the **effective force on the single DOF** — the caller has already lumped
/// the loaded tributary area into it (force = pressure × area), per the
/// reused [`solve_sdof_blast`] convention.
///
/// # Errors
///
/// - [`SurvivabilityError::InvalidParameter`] for a non-positive `mass` or
///   `stiffness`, or a non-positive `yield_deflection_m` when supplied.
/// - [`SurvivabilityError::Transient`] if the underlying [`valenx_fem`] solve
///   fails.
pub fn sdof_response(
    mass: f64,
    damping: f64,
    stiffness: f64,
    yield_deflection_m: Option<f64>,
    pulse: &FriedlanderPulse,
    controls: &TransientControls,
) -> Result<SdofResponse, SurvivabilityError> {
    let m = require_pos("effective mass m", mass)?;
    let k = require_pos("effective stiffness k", stiffness)?;
    if !(damping.is_finite() && damping >= 0.0) {
        return Err(SurvivabilityError::InvalidParameter(format!(
            "damping c must be finite and >= 0, got {damping}"
        )));
    }
    if let Some(xy) = yield_deflection_m {
        require_pos("yield deflection", xy)?;
    }

    let resp = solve_sdof_blast(m, damping, k, pulse, controls)?;
    let peak = resp.peak_abs_displacement(0);
    let peak_acc = resp.peak_abs_acceleration(0);

    // k > 0 (checked) ⇒ static deflection divide is safe.
    let x_st = pulse.peak_overpressure / k;
    let daf = if x_st > 0.0 { peak / x_st } else { 0.0 };

    let ductility = yield_deflection_m.map(|xy| peak / xy);

    Ok(SdofResponse {
        peak_deflection_m: peak,
        peak_acceleration_m_s2: peak_acc,
        ductility_ratio: ductility,
        static_deflection_m: x_st,
        dynamic_amplification: daf,
    })
}

/// One point on a pressure–impulse iso-damage curve: the threshold peak
/// pressure and specific impulse that *together* just reach the damage limit.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PiPoint {
    /// Threshold peak pressure `P` (Pa) at this point on the curve.
    pub pressure_pa: f64,
    /// Threshold specific impulse `i` (Pa·s) at this point on the curve.
    pub impulse_pa_s: f64,
}

/// A **pressure–impulse (P–I) diagram** for an elastic SDOF protective element
/// at a chosen deflection (damage) limit.
///
/// ## Model
///
/// For an undamped SDOF of effective mass `m`, stiffness `k`, natural frequency
/// `ω = √(k/m)`, and an elastic deflection limit `x_lim`, energy balance gives
/// the two asymptotes exactly:
///
/// - **Impulsive asymptote** (load far shorter than the period — all the
///   impulse is delivered before the structure moves): kinetic energy
///   `½ m v² = ½ m (i/m)²` equals strain energy `½ k x_lim²`, so the **minimum
///   damaging impulse** is
///
///   ```text
///     i_min = x_lim · √(k · m) = m · ω · x_lim
///   ```
///
///   a **vertical** line independent of pressure.
///
/// - **Quasi-static asymptote** (load far longer than the period — the peak
///   pressure is effectively static): the structure reaches the limit when the
///   *static* deflection `P·A/k` (per unit effective area, here normalized so
///   the resistance `R = k·x_lim`) equals `x_lim`, i.e. the **minimum damaging
///   pressure** is
///
///   ```text
///     P_min = k · x_lim          (≡ the element's static resistance R)
///   ```
///
///   a **horizontal** line independent of impulse.
///
/// The full iso-damage curve interpolates between them with a rectangular
/// hyperbola asymptotic to both lines (Baker et al., *Explosion Hazards and
/// Evaluation*, Elsevier 1983; Smith & Hetherington, *Blast and Ballistic
/// Loading of Structures*, 1994):
///
/// ```text
///   (P − P_min)(i − i_min) = ½ · P_min · i_min          (P > P_min, i > i_min)
/// ```
///
/// As `i → ∞` the left factor must vanish, so `P → P_min` (the quasi-static
/// asymptote); as `P → ∞`, `i → i_min` (the impulsive asymptote). The constant
/// `½·P_min·i_min` fixes the curvature of the knee between them.
/// Research/educational, validation-pending.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct PressureImpulseDiagram {
    /// Quasi-static (horizontal) asymptote — minimum damaging pressure
    /// `P_min = k·x_lim` (Pa).
    pub pressure_asymptote_pa: f64,
    /// Impulsive (vertical) asymptote — minimum damaging impulse
    /// `i_min = x_lim·√(k·m)` (Pa·s).
    pub impulse_asymptote_pa_s: f64,
}

impl PressureImpulseDiagram {
    /// Build the P–I diagram for an elastic SDOF element of effective `mass`
    /// (kg) and `stiffness` (N/m) at deflection limit `deflection_limit_m` (m).
    ///
    /// # Errors
    ///
    /// [`SurvivabilityError::InvalidParameter`] if any input is not
    /// finite-and-positive.
    pub fn elastic(
        mass: f64,
        stiffness: f64,
        deflection_limit_m: f64,
    ) -> Result<PressureImpulseDiagram, SurvivabilityError> {
        let m = require_pos("effective mass m", mass)?;
        let k = require_pos("effective stiffness k", stiffness)?;
        let x_lim = require_pos("deflection limit", deflection_limit_m)?;

        // Both asymptotes are products/roots of positive quantities — safe.
        let p_min = k * x_lim;
        let i_min = x_lim * (k * m).sqrt();

        Ok(PressureImpulseDiagram {
            pressure_asymptote_pa: p_min,
            impulse_asymptote_pa_s: i_min,
        })
    }

    /// The threshold impulse on the iso-damage curve at a given peak `pressure`
    /// (Pa). For `pressure > P_min` the hyperbola gives a finite impulse that
    /// **decreases toward `i_min`** as pressure grows. Returns `None` if
    /// `pressure ≤ P_min` (at or left of the quasi-static asymptote — no finite
    /// damaging impulse exists there; the point is safe on pressure alone).
    pub fn impulse_at_pressure(&self, pressure_pa: f64) -> Option<f64> {
        let p = pressure_pa;
        let p_min = self.pressure_asymptote_pa;
        let i_min = self.impulse_asymptote_pa_s;
        let denom = p - p_min;
        if !p.is_finite() || denom <= 0.0 {
            return None;
        }
        // (P − P_min)(i − i_min) = ½ P_min i_min  ⇒  i = i_min + ½P_min·i_min/denom
        Some(i_min + 0.5 * p_min * i_min / denom)
    }

    /// The threshold pressure on the iso-damage curve at a given specific
    /// `impulse` (Pa·s). Symmetric to [`impulse_at_pressure`](Self::impulse_at_pressure).
    /// `None` if `impulse ≤ i_min` (at or below the impulsive asymptote).
    pub fn pressure_at_impulse(&self, impulse_pa_s: f64) -> Option<f64> {
        let i = impulse_pa_s;
        let p_min = self.pressure_asymptote_pa;
        let i_min = self.impulse_asymptote_pa_s;
        let denom = i - i_min;
        if !i.is_finite() || denom <= 0.0 {
            return None;
        }
        Some(p_min + 0.5 * p_min * i_min / denom)
    }

    /// Is a load `(pressure, impulse)` **safe** (at or below the iso-damage
    /// curve)? A point is damaging when it lies above-and-right of the curve,
    /// i.e. when `(P − P_min)(i − i_min) ≥ ½ P_min i_min` with both factors
    /// positive. Anything with `P ≤ P_min` or `i ≤ i_min` is below an asymptote
    /// and therefore safe.
    pub fn is_safe(&self, pressure_pa: f64, impulse_pa_s: f64) -> bool {
        let p_min = self.pressure_asymptote_pa;
        let i_min = self.impulse_asymptote_pa_s;
        let fp = pressure_pa - p_min;
        let fi = impulse_pa_s - i_min;
        if fp <= 0.0 || fi <= 0.0 {
            return true;
        }
        fp * fi < 0.5 * p_min * i_min
    }

    /// Sample `n` points along the iso-damage curve, logarithmically spaced in
    /// pressure from just above the quasi-static asymptote up to `factor`×P_min,
    /// for plotting. Returns at least the two endpoints. `n` is clamped to
    /// `[2, 4096]`; `factor` to `≥ 1.0001`.
    pub fn curve(&self, n: usize, factor: f64) -> Vec<PiPoint> {
        let n = n.clamp(2, 4096);
        let factor = if factor.is_finite() && factor > 1.0001 {
            factor
        } else {
            100.0
        };
        let p_min = self.pressure_asymptote_pa;
        // Start a hair above the asymptote so the impulse is finite.
        let p_start = p_min * 1.01;
        let p_end = p_min * factor;
        let ln0 = p_start.ln();
        let ln1 = p_end.ln();
        let mut out = Vec::with_capacity(n);
        for j in 0..n {
            let frac = j as f64 / (n as f64 - 1.0);
            let p = (ln0 + frac * (ln1 - ln0)).exp();
            if let Some(i) = self.impulse_at_pressure(p) {
                out.push(PiPoint {
                    pressure_pa: p,
                    impulse_pa_s: i,
                });
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blast::BlastLoad;

    #[test]
    fn undamped_impulse_response_peak_matches_analytic() {
        // PIN (pure SDOF impulse response): an undamped SDOF given an impulse I
        // over a time short vs its period starts with velocity v = I/m and then
        // oscillates as x(t) = (I/(mω))·sin(ωt), so the peak deflection is
        //     x_max = (I/m)/ω = v/ω.
        // We drive the reused valenx-fem Newmark integrator directly with a
        // short *rectangular* force pulse (constant F over [0, t_d], zero after)
        // so the delivered impulse is exactly F·t_d with no suction tail — this
        // isolates the SDOF impulse-response physics. (The Friedlander pulse is
        // exercised end-to-end in `end_to_end_blast_to_response`; its full form
        // carries a negative phase and so is not a pure impulse.)
        use nalgebra::{DMatrix, DVector};

        let m: f64 = 10.0;
        let k: f64 = 4.0e6; // ω = sqrt(k/m) = sqrt(4e5) ≈ 632.456 rad/s
        let omega = (k / m).sqrt();
        let period = 2.0 * std::f64::consts::PI / omega; // ≈ 9.93 ms

        let f_peak = 1.0e6;
        let t_d = period / 500.0; // genuinely impulsive
        let impulse = f_peak * t_d; // exact rectangular impulse

        let dt = t_d / 20.0;
        let n_steps = (0.3 * period / dt) as usize; // past the T/4 peak
        let controls = TransientControls {
            dt,
            n_steps,
            newmark: valenx_fem::NewmarkBeta::average_acceleration(),
        };

        let mass = DMatrix::from_element(1, 1, m);
        let damp = DMatrix::from_element(1, 1, 0.0); // undamped
        let stiff = DMatrix::from_element(1, 1, k);
        let u0 = DVector::zeros(1);
        let v0 = DVector::zeros(1);
        let resp = valenx_fem::solve_transient_response(
            &mass,
            &damp,
            &stiff,
            &u0,
            &v0,
            |t| DVector::from_element(1, if (0.0..t_d).contains(&t) { f_peak } else { 0.0 }),
            &controls,
        )
        .unwrap();

        let peak = resp.peak_abs_displacement(0);
        let analytic_peak = (impulse / m) / omega;
        let rel = (peak - analytic_peak).abs() / analytic_peak;
        assert!(
            rel < 0.05,
            "impulse-response peak off: got {peak}, analytic {analytic_peak}, rel {rel}"
        );
    }

    #[test]
    fn pi_diagram_asymptotes() {
        // PIN: the two asymptotes are i_min = x_lim·√(km) and P_min = k·x_lim.
        let m = 50.0;
        let k = 2.0e6;
        let x_lim = 0.01;
        let pi = PressureImpulseDiagram::elastic(m, k, x_lim).unwrap();

        let expect_p = k * x_lim; // 20 000 Pa
        let expect_i = x_lim * (k * m).sqrt(); // 0.01·sqrt(1e8)=0.01·1e4=100 Pa·s
        assert!((pi.pressure_asymptote_pa - expect_p).abs() < 1e-6);
        assert!((pi.impulse_asymptote_pa_s - expect_i).abs() < 1e-6);
    }

    #[test]
    fn pi_curve_approaches_asymptotes() {
        let pi = PressureImpulseDiagram::elastic(50.0, 2.0e6, 0.01).unwrap();
        let p_min = pi.pressure_asymptote_pa;
        let i_min = pi.impulse_asymptote_pa_s;

        // Quasi-static limit: as impulse → ∞, threshold pressure → P_min⁺.
        let p_hi_i = pi.pressure_at_impulse(i_min * 1.0e6).unwrap();
        assert!(
            (p_hi_i - p_min).abs() / p_min < 1e-3,
            "pressure should approach P_min for huge impulse: {p_hi_i} vs {p_min}"
        );
        // Impulsive limit: as pressure → ∞, threshold impulse → i_min⁺.
        let i_hi_p = pi.impulse_at_pressure(p_min * 1.0e6).unwrap();
        assert!(
            (i_hi_p - i_min).abs() / i_min < 1e-3,
            "impulse should approach i_min for huge pressure: {i_hi_p} vs {i_min}"
        );
    }

    #[test]
    fn pi_curve_monotone_and_safe_region() {
        let pi = PressureImpulseDiagram::elastic(50.0, 2.0e6, 0.01).unwrap();
        // A point well inside both half-asymptotes is safe.
        assert!(pi.is_safe(
            pi.pressure_asymptote_pa * 0.4,
            pi.impulse_asymptote_pa_s * 0.4
        ));
        // A point far above-and-right of the knee is unsafe (damaging).
        assert!(!pi.is_safe(
            pi.pressure_asymptote_pa * 10.0,
            pi.impulse_asymptote_pa_s * 10.0
        ));
        // The plotting curve is non-empty and impulse decreases with pressure.
        let curve = pi.curve(50, 100.0);
        assert!(curve.len() >= 2);
        for w in curve.windows(2) {
            assert!(w[1].pressure_pa > w[0].pressure_pa);
            assert!(w[1].impulse_pa_s <= w[0].impulse_pa_s + 1e-9);
        }
    }

    #[test]
    fn response_degenerate_inputs_error() {
        let pulse = FriedlanderPulse::new(1.0e5, 1.0e-3, 1.0).unwrap();
        let c = TransientControls::default();
        assert!(sdof_response(0.0, 0.0, 1.0e6, None, &pulse, &c).is_err()); // zero mass
        assert!(sdof_response(10.0, 0.0, 0.0, None, &pulse, &c).is_err()); // zero stiff
        assert!(sdof_response(10.0, -1.0, 1.0e6, None, &pulse, &c).is_err()); // neg damping
        assert!(PressureImpulseDiagram::elastic(0.0, 1.0, 1.0).is_err());
        assert!(PressureImpulseDiagram::elastic(1.0, 0.0, 1.0).is_err());
        assert!(PressureImpulseDiagram::elastic(1.0, 1.0, 0.0).is_err());
    }

    #[test]
    fn end_to_end_blast_to_response() {
        // A realistic protective screen: 100 kg TNT at 20 m onto a wall panel
        // reduced to an effective SDOF. (Force already lumped into peak for the
        // SDOF convention; here we just confirm the pipeline runs and is sane.)
        let load = BlastLoad::tnt_free_air(100.0, 20.0).unwrap();
        let pulse = load.friedlander().unwrap();
        let c = TransientControls {
            dt: load.positive_duration_s / 50.0,
            n_steps: 4000,
            newmark: valenx_fem::NewmarkBeta::average_acceleration(),
        };
        let resp = sdof_response(500.0, 100.0, 5.0e6, Some(0.02), &pulse, &c).unwrap();
        assert!(resp.peak_deflection_m > 0.0);
        assert!(resp.dynamic_amplification > 0.0);
        assert!(resp.ductility_ratio.unwrap() > 0.0);

        let json = serde_json::to_string(&resp).unwrap();
        let _back: SdofResponse = serde_json::from_str(&json).unwrap();
    }
}
