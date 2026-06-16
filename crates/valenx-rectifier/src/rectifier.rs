//! Closed-form averages, RMS, and ripple factors for ideal diode
//! rectifiers driven by a sinusoidal mains.
//!
//! All quantities are referred to the peak of the (ideal, lossless)
//! sinusoidal input `v(t) = Vpeak * sin(omega t)`. Diodes are treated as
//! ideal switches (zero forward drop, zero reverse leakage), so the
//! results are the textbook first-order figures, not measurements of a
//! real bridge.
//!
//! # Definitions
//!
//! For a waveform `v(t)` over one period `T`:
//!
//! - average (DC) value: `Vdc = (1/T) * integral over T of v(t) dt`,
//! - root-mean-square value: `Vrms = sqrt((1/T) * integral over T of v(t)^2 dt)`,
//! - ripple factor: `r = Vac_rms / Vdc = sqrt((Vrms/Vdc)^2 - 1)`,
//!   the ratio of the RMS of the AC (ripple) component to the DC value.

use crate::error::RectifierError;

/// Which diode topology is being analysed.
///
/// The two variants differ by a factor of two in both average output and
/// (for a capacitor filter) effective ripple frequency: a full-wave
/// rectifier conducts on both half-cycles of the mains.
#[derive(Copy, Clone, Debug, Eq, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Topology {
    /// Single diode; conducts on one half-cycle per mains period.
    HalfWave,
    /// Bridge / centre-tapped; conducts on both half-cycles per period.
    FullWave,
}

impl Topology {
    /// Output-ripple frequency as a multiple of the mains frequency.
    ///
    /// A half-wave rectifier produces one output hump per mains cycle, so
    /// its ripple is at the mains frequency (multiplier `1`). A full-wave
    /// rectifier produces two humps per cycle, doubling the ripple
    /// frequency (multiplier `2`). This is the multiplier that turns the
    /// mains frequency into the `f` used in [`capacitor_ripple_pp`].
    pub fn ripple_frequency_multiplier(self) -> f64 {
        match self {
            Topology::HalfWave => 1.0,
            Topology::FullWave => 2.0,
        }
    }
}

/// Average (DC) output voltage of an ideal **half-wave** rectifier.
///
/// `Vdc = Vpeak / pi`.
///
/// Over one period the single diode passes the positive half-sine and
/// blocks the negative half, so the mean is `Vpeak/pi ~= 0.3183 Vpeak`.
///
/// # Errors
///
/// Returns [`RectifierError`] if `v_peak` is not finite and positive.
pub fn half_wave_vdc(v_peak: f64) -> Result<f64, RectifierError> {
    let v_peak = RectifierError::positive("v_peak", v_peak)?;
    Ok(v_peak / std::f64::consts::PI)
}

/// Average (DC) output voltage of an ideal **full-wave** rectifier.
///
/// `Vdc = 2 * Vpeak / pi`.
///
/// Both half-cycles are passed (the negative one inverted), so the mean
/// is exactly twice the half-wave value, `2 Vpeak/pi ~= 0.6366 Vpeak`.
///
/// # Errors
///
/// Returns [`RectifierError`] if `v_peak` is not finite and positive.
pub fn full_wave_vdc(v_peak: f64) -> Result<f64, RectifierError> {
    let v_peak = RectifierError::positive("v_peak", v_peak)?;
    Ok(2.0 * v_peak / std::f64::consts::PI)
}

/// RMS output voltage of an ideal **half-wave** rectifier.
///
/// `Vrms = Vpeak / 2`.
///
/// The squared waveform is non-zero only over the conducting half-cycle,
/// halving the mean-square relative to a full sine, hence `Vpeak/2`.
///
/// # Errors
///
/// Returns [`RectifierError`] if `v_peak` is not finite and positive.
pub fn half_wave_vrms(v_peak: f64) -> Result<f64, RectifierError> {
    let v_peak = RectifierError::positive("v_peak", v_peak)?;
    Ok(v_peak / 2.0)
}

/// RMS output voltage of an ideal **full-wave** rectifier.
///
/// `Vrms = Vpeak / sqrt(2)`.
///
/// Rectification does not change the mean-square (it only flips sign on
/// alternate half-cycles), so the RMS equals that of the underlying sine,
/// `Vpeak/sqrt(2) ~= 0.7071 Vpeak`.
///
/// # Errors
///
/// Returns [`RectifierError`] if `v_peak` is not finite and positive.
pub fn full_wave_vrms(v_peak: f64) -> Result<f64, RectifierError> {
    let v_peak = RectifierError::positive("v_peak", v_peak)?;
    Ok(v_peak / std::f64::consts::SQRT_2)
}

/// Average (DC) output for the given [`Topology`].
///
/// Dispatches to [`half_wave_vdc`] or [`full_wave_vdc`].
///
/// # Errors
///
/// Returns [`RectifierError`] if `v_peak` is not finite and positive.
pub fn vdc(topology: Topology, v_peak: f64) -> Result<f64, RectifierError> {
    match topology {
        Topology::HalfWave => half_wave_vdc(v_peak),
        Topology::FullWave => full_wave_vdc(v_peak),
    }
}

/// RMS output for the given [`Topology`].
///
/// Dispatches to [`half_wave_vrms`] or [`full_wave_vrms`].
///
/// # Errors
///
/// Returns [`RectifierError`] if `v_peak` is not finite and positive.
pub fn vrms(topology: Topology, v_peak: f64) -> Result<f64, RectifierError> {
    match topology {
        Topology::HalfWave => half_wave_vrms(v_peak),
        Topology::FullWave => full_wave_vrms(v_peak),
    }
}

/// Ripple factor from an RMS / DC pair.
///
/// `r = sqrt((Vrms / Vdc)^2 - 1)`.
///
/// This is the ratio of the RMS of the AC (ripple) component to the DC
/// component; it is independent of `Vpeak` for an ideal rectifier and so
/// is a pure dimensionless figure of merit (smaller is smoother).
///
/// # Errors
///
/// Returns [`RectifierError`] if either argument is not finite, if
/// `v_dc` is not strictly positive, or if `v_rms` is negative. The
/// physical relation guarantees `Vrms >= Vdc`; if a caller passes a pair
/// with `Vrms < Vdc` the radicand would be negative, so the result is
/// clamped to `0.0` (a perfectly smooth output) rather than yielding
/// `NaN`.
pub fn ripple_factor_from(v_rms: f64, v_dc: f64) -> Result<f64, RectifierError> {
    let v_rms = RectifierError::non_negative("v_rms", v_rms)?;
    let v_dc = RectifierError::positive("v_dc", v_dc)?;
    let ratio = v_rms / v_dc;
    let radicand = ratio * ratio - 1.0;
    Ok(if radicand <= 0.0 {
        0.0
    } else {
        radicand.sqrt()
    })
}

/// Ripple factor of an ideal rectifier of the given [`Topology`].
///
/// Convenience wrapper that computes [`vrms`] and [`vdc`] for a probe
/// `Vpeak` and feeds them to [`ripple_factor_from`]. Because the factor
/// is independent of `Vpeak`, the constant `1.0` is used internally.
///
/// The closed-form values are `r_half = sqrt(pi^2/4 - 1) ~= 1.211` and
/// `r_full = sqrt(pi^2/8 - 1) ~= 0.483`, so a half-wave output is much
/// rougher than a full-wave one.
///
/// # Errors
///
/// Never fails for the internal probe value, but returns the
/// [`RectifierError`] type for signature symmetry with the rest of the
/// module.
pub fn ripple_factor(topology: Topology) -> Result<f64, RectifierError> {
    let v_peak = 1.0;
    let v_rms = vrms(topology, v_peak)?;
    let v_dc = vdc(topology, v_peak)?;
    ripple_factor_from(v_rms, v_dc)
}

/// Form factor from an RMS / DC pair.
///
/// `FF = Vrms / Vdc`.
///
/// The form factor measures how far a rectified waveform departs from
/// pure DC: it is `1` for ideal smooth DC and grows as the output gets
/// rougher. It is the companion ratio to the [`ripple_factor`], tied to
/// it exactly by `r = sqrt(FF^2 - 1)` and `FF = sqrt(r^2 + 1)`, and is
/// independent of `Vpeak` for an ideal rectifier.
///
/// # Errors
///
/// Returns [`RectifierError`] if either argument is not finite, if `v_dc`
/// is not strictly positive, or if `v_rms` is negative.
pub fn form_factor_from(v_rms: f64, v_dc: f64) -> Result<f64, RectifierError> {
    let v_rms = RectifierError::non_negative("v_rms", v_rms)?;
    let v_dc = RectifierError::positive("v_dc", v_dc)?;
    Ok(v_rms / v_dc)
}

/// Form factor of an ideal rectifier of the given [`Topology`].
///
/// Convenience wrapper that computes [`vrms`] and [`vdc`] for a probe
/// `Vpeak` and feeds them to [`form_factor_from`]. Because the factor is
/// independent of `Vpeak`, the constant `1.0` is used internally.
///
/// The closed-form values are `FF_half = pi/2 ~= 1.571` and
/// `FF_full = pi/(2*sqrt(2)) ~= 1.111`, so a full-wave output is far
/// closer to smooth DC than a half-wave one.
///
/// # Errors
///
/// Returns the [`RectifierError`] type for signature symmetry with the
/// rest of the module; never fails for the internal probe value.
pub fn form_factor(topology: Topology) -> Result<f64, RectifierError> {
    let v_peak = 1.0;
    let v_rms = vrms(topology, v_peak)?;
    let v_dc = vdc(topology, v_peak)?;
    form_factor_from(v_rms, v_dc)
}

/// Rectification efficiency from an RMS / DC pair.
///
/// `eta = (Vdc / Vrms)^2`.
///
/// The ratio of the DC output power `Vdc^2 / R` to the total power
/// `Vrms^2 / R` delivered to a resistive load — the fraction of the load
/// power that is useful DC. It is the reciprocal of the squared form
/// factor (`eta = 1 / FF^2 = 1 / (1 + r^2)`), independent of `Vpeak`, and
/// lies in `(0, 1]`.
///
/// # Errors
///
/// Returns [`RectifierError`] if either argument is not finite, if
/// `v_rms` is not strictly positive, or if `v_dc` is negative.
pub fn rectification_efficiency_from(v_rms: f64, v_dc: f64) -> Result<f64, RectifierError> {
    let v_rms = RectifierError::positive("v_rms", v_rms)?;
    let v_dc = RectifierError::non_negative("v_dc", v_dc)?;
    let ratio = v_dc / v_rms;
    Ok(ratio * ratio)
}

/// Rectification efficiency of an ideal rectifier of the given
/// [`Topology`].
///
/// Convenience wrapper feeding [`vrms`] and [`vdc`] for a probe `Vpeak`
/// to [`rectification_efficiency_from`]. The famous closed-form values
/// are `eta_half = 4/pi^2 ~= 0.406` (40.6%) and
/// `eta_full = 8/pi^2 ~= 0.811` (81.1%): a full-wave rectifier delivers
/// twice the useful DC-power fraction of a half-wave one.
///
/// # Errors
///
/// Returns the [`RectifierError`] type for signature symmetry; never
/// fails for the internal probe value.
pub fn rectification_efficiency(topology: Topology) -> Result<f64, RectifierError> {
    let v_peak = 1.0;
    let v_rms = vrms(topology, v_peak)?;
    let v_dc = vdc(topology, v_peak)?;
    rectification_efficiency_from(v_rms, v_dc)
}

/// Peak-to-peak ripple voltage of a capacitor-input filter.
///
/// `Vr = I / (f * C)`.
///
/// First-order approximation in which the load draws a constant current
/// `I` (amperes) while the reservoir capacitor `C` (farads) discharges
/// linearly between conduction pulses arriving at frequency `f` (hertz).
/// For a full-wave rectifier on a mains of frequency `f_mains`, use
/// `f = 2 * f_mains` (see [`Topology::ripple_frequency_multiplier`]); for
/// a half-wave rectifier use `f = f_mains`. Larger `C` (or higher `f`)
/// yields proportionally smaller ripple.
///
/// # Errors
///
/// Returns [`RectifierError`] if `load_current_a` is negative or any
/// argument is non-finite, or if `freq_hz` or `cap_farads` is not
/// strictly positive.
pub fn capacitor_ripple_pp(
    load_current_a: f64,
    freq_hz: f64,
    cap_farads: f64,
) -> Result<f64, RectifierError> {
    let i = RectifierError::non_negative("load_current_a", load_current_a)?;
    let f = RectifierError::positive("freq_hz", freq_hz)?;
    let c = RectifierError::positive("cap_farads", cap_farads)?;
    Ok(i / (f * c))
}

/// Peak-to-peak capacitor ripple given a [`Topology`] and the **mains**
/// frequency.
///
/// Computes the effective ripple frequency
/// `f = mains_freq_hz * topology.ripple_frequency_multiplier()` and
/// forwards to [`capacitor_ripple_pp`]. This is the convenient entry
/// point when you have the line frequency (e.g. 50 or 60 Hz) rather than
/// the already-doubled ripple frequency.
///
/// # Errors
///
/// Propagates the validation errors of [`capacitor_ripple_pp`]; also
/// returns [`RectifierError`] if `mains_freq_hz` is not finite and
/// positive.
pub fn capacitor_ripple_pp_for(
    topology: Topology,
    load_current_a: f64,
    mains_freq_hz: f64,
    cap_farads: f64,
) -> Result<f64, RectifierError> {
    let mains = RectifierError::positive("mains_freq_hz", mains_freq_hz)?;
    let f = mains * topology.ripple_frequency_multiplier();
    capacitor_ripple_pp(load_current_a, f, cap_farads)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const EPS: f64 = 1e-12;

    // ---- Average (DC) values --------------------------------------------

    #[test]
    fn half_wave_dc_is_vpeak_over_pi() {
        let v = half_wave_vdc(10.0).expect("valid peak");
        assert!((v - 10.0 / PI).abs() < EPS, "got {v}");
        // Numeric sanity: ~3.1831 for a 10 V peak.
        assert!((v - 3.183_098_861_837_907).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn full_wave_dc_is_two_vpeak_over_pi() {
        let v = full_wave_vdc(10.0).expect("valid peak");
        assert!((v - 2.0 * 10.0 / PI).abs() < EPS, "got {v}");
    }

    #[test]
    fn full_wave_dc_is_exactly_twice_half_wave_dc() {
        let peak = 17.3;
        let half = half_wave_vdc(peak).expect("valid");
        let full = full_wave_vdc(peak).expect("valid");
        assert!((full - 2.0 * half).abs() < EPS, "half={half} full={full}");
    }

    // ---- RMS values ------------------------------------------------------

    #[test]
    fn half_wave_rms_is_vpeak_over_two() {
        let v = half_wave_vrms(10.0).expect("valid");
        assert!((v - 5.0).abs() < EPS, "got {v}");
    }

    #[test]
    fn full_wave_rms_is_vpeak_over_sqrt2() {
        let v = full_wave_vrms(10.0).expect("valid");
        assert!((v - 10.0 / std::f64::consts::SQRT_2).abs() < EPS, "got {v}");
        assert!((v - 7.071_067_811_865_476).abs() < 1e-9, "got {v}");
    }

    #[test]
    fn full_wave_rms_exceeds_half_wave_rms() {
        let peak = 230.0_f64 * std::f64::consts::SQRT_2; // 230 V RMS mains peak
        let half = half_wave_vrms(peak).expect("valid");
        let full = full_wave_vrms(peak).expect("valid");
        assert!(full > half, "half={half} full={full}");
    }

    // ---- Dispatchers agree with the explicit functions ------------------

    #[test]
    fn dispatchers_match_explicit() {
        let peak = 12.0;
        assert!(
            (vdc(Topology::HalfWave, peak).unwrap() - half_wave_vdc(peak).unwrap()).abs() < EPS
        );
        assert!(
            (vdc(Topology::FullWave, peak).unwrap() - full_wave_vdc(peak).unwrap()).abs() < EPS
        );
        assert!(
            (vrms(Topology::HalfWave, peak).unwrap() - half_wave_vrms(peak).unwrap()).abs() < EPS
        );
        assert!(
            (vrms(Topology::FullWave, peak).unwrap() - full_wave_vrms(peak).unwrap()).abs() < EPS
        );
    }

    // ---- Ripple factor ---------------------------------------------------

    #[test]
    fn ripple_factor_half_wave_is_about_1_21() {
        // Closed form: sqrt((Vrms/Vdc)^2 - 1) with Vrms=Vpeak/2,
        // Vdc=Vpeak/pi  ->  sqrt(pi^2/4 - 1).
        let r = ripple_factor(Topology::HalfWave).expect("valid");
        let expected = (PI * PI / 4.0 - 1.0).sqrt();
        assert!((r - expected).abs() < EPS, "got {r}");
        // Standard textbook value ~1.21.
        assert!((r - 1.211_363_322_984_619).abs() < 1e-9, "got {r}");
    }

    #[test]
    fn ripple_factor_full_wave_is_about_0_48() {
        // sqrt((Vpeak/sqrt2 / (2 Vpeak/pi))^2 - 1) = sqrt(pi^2/8 - 1).
        let r = ripple_factor(Topology::FullWave).expect("valid");
        let expected = (PI * PI / 8.0 - 1.0).sqrt();
        assert!((r - expected).abs() < EPS, "got {r}");
        // Standard textbook value ~0.48.
        assert!((r - 0.483_425_847_608_678).abs() < 1e-9, "got {r}");
    }

    #[test]
    fn ripple_factor_half_wave_exceeds_full_wave() {
        let half = ripple_factor(Topology::HalfWave).expect("valid");
        let full = ripple_factor(Topology::FullWave).expect("valid");
        assert!(half > full, "half={half} full={full}");
    }

    #[test]
    fn ripple_factor_is_independent_of_peak() {
        // Compute from explicit pairs at two very different peaks; the
        // dimensionless factor must be identical.
        let r_small =
            ripple_factor_from(half_wave_vrms(1.0).unwrap(), half_wave_vdc(1.0).unwrap()).unwrap();
        let r_big = ripple_factor_from(
            half_wave_vrms(1.0e6).unwrap(),
            half_wave_vdc(1.0e6).unwrap(),
        )
        .unwrap();
        assert!((r_small - r_big).abs() < EPS, "small={r_small} big={r_big}");
    }

    #[test]
    fn ripple_factor_clamps_when_rms_below_dc() {
        // Non-physical pair (Vrms < Vdc) clamps to a smooth 0.0 instead
        // of returning NaN from a negative radicand.
        let r = ripple_factor_from(1.0, 2.0).expect("valid inputs");
        assert!(r.abs() < EPS, "got {r}");
    }

    // ---- Form factor and rectification efficiency -----------------------

    #[test]
    fn form_factor_half_wave_is_pi_over_two() {
        let ff = form_factor(Topology::HalfWave).expect("valid");
        // Closed form FF_half = pi/2 ~= 1.5708 (the standard value).
        assert!((ff - PI / 2.0).abs() < EPS, "got {ff}");
        assert!((ff - std::f64::consts::FRAC_PI_2).abs() < EPS, "got {ff}");
    }

    #[test]
    fn form_factor_full_wave_is_pi_over_two_root_two() {
        let ff = form_factor(Topology::FullWave).expect("valid");
        let expected = PI / (2.0 * std::f64::consts::SQRT_2);
        assert!((ff - expected).abs() < EPS, "got {ff}");
        assert!((ff - 1.110_720_734_539_591).abs() < 1e-9, "got {ff}");
    }

    #[test]
    fn rectification_efficiency_half_wave_is_about_40_percent() {
        let eta = rectification_efficiency(Topology::HalfWave).expect("valid");
        assert!((eta - 4.0 / (PI * PI)).abs() < EPS, "got {eta}");
        // Standard textbook value ~0.406 (40.6%).
        assert!((eta - 0.405_284_734_569_351).abs() < 1e-9, "got {eta}");
    }

    #[test]
    fn rectification_efficiency_full_wave_is_about_81_percent() {
        let eta = rectification_efficiency(Topology::FullWave).expect("valid");
        assert!((eta - 8.0 / (PI * PI)).abs() < EPS, "got {eta}");
        // Standard textbook value ~0.811 (81.1%), exactly twice the
        // half-wave figure.
        assert!((eta - 0.810_569_469_138_702).abs() < 1e-9, "got {eta}");
    }

    #[test]
    fn form_factor_and_ripple_factor_satisfy_the_identity() {
        // r = sqrt(FF^2 - 1), i.e. FF = sqrt(r^2 + 1), for both topologies.
        for t in [Topology::HalfWave, Topology::FullWave] {
            let ff = form_factor(t).unwrap();
            let r = ripple_factor(t).unwrap();
            assert!((r - (ff * ff - 1.0).sqrt()).abs() < EPS, "{t:?}: r vs FF");
            assert!(((r * r + 1.0).sqrt() - ff).abs() < EPS, "{t:?}: FF vs r");
        }
    }

    #[test]
    fn efficiency_is_reciprocal_of_squared_form_factor() {
        // eta = 1 / FF^2 = 1 / (1 + r^2).
        for t in [Topology::HalfWave, Topology::FullWave] {
            let eta = rectification_efficiency(t).unwrap();
            let ff = form_factor(t).unwrap();
            let r = ripple_factor(t).unwrap();
            assert!((eta - 1.0 / (ff * ff)).abs() < EPS, "{t:?}: eta vs 1/FF^2");
            assert!(
                (eta - 1.0 / (1.0 + r * r)).abs() < EPS,
                "{t:?}: eta vs 1/(1+r^2)"
            );
        }
    }

    #[test]
    fn full_wave_is_smoother_and_more_efficient_than_half_wave() {
        // Lower form factor (closer to 1) and higher efficiency.
        assert!(
            form_factor(Topology::FullWave).unwrap() < form_factor(Topology::HalfWave).unwrap()
        );
        assert!(
            rectification_efficiency(Topology::FullWave).unwrap()
                > rectification_efficiency(Topology::HalfWave).unwrap()
        );
    }

    #[test]
    fn efficiency_lies_in_zero_to_one() {
        for t in [Topology::HalfWave, Topology::FullWave] {
            let eta = rectification_efficiency(t).unwrap();
            assert!(eta > 0.0 && eta <= 1.0, "{t:?}: eta = {eta}");
        }
    }

    #[test]
    fn form_factor_is_independent_of_peak() {
        let small =
            form_factor_from(half_wave_vrms(1.0).unwrap(), half_wave_vdc(1.0).unwrap()).unwrap();
        let big = form_factor_from(
            half_wave_vrms(1.0e6).unwrap(),
            half_wave_vdc(1.0e6).unwrap(),
        )
        .unwrap();
        assert!((small - big).abs() < EPS, "small={small} big={big}");
    }

    #[test]
    fn form_factor_and_efficiency_reject_bad_inputs() {
        assert!(form_factor_from(1.0, 0.0).is_err()); // zero dc (denominator)
        assert!(form_factor_from(-1.0, 1.0).is_err()); // negative rms
        assert!(rectification_efficiency_from(0.0, 1.0).is_err()); // zero rms (denominator)
        assert!(rectification_efficiency_from(1.0, -1.0).is_err()); // negative dc
        assert!(rectification_efficiency_from(f64::NAN, 1.0).is_err()); // non-finite
    }

    // ---- Capacitor-filter ripple ----------------------------------------

    #[test]
    fn capacitor_ripple_is_i_over_fc() {
        // 1 A load, 100 Hz (full-wave on 50 Hz mains), 1000 uF.
        let vr = capacitor_ripple_pp(1.0, 100.0, 1.0e-3).expect("valid");
        assert!((vr - 10.0).abs() < EPS, "got {vr}");
    }

    #[test]
    fn larger_capacitor_gives_less_ripple() {
        let small = capacitor_ripple_pp(0.5, 120.0, 470.0e-6).expect("valid");
        let large = capacitor_ripple_pp(0.5, 120.0, 4700.0e-6).expect("valid");
        assert!(large < small, "small_c={small} large_c={large}");
        // 10x the capacitance -> 1/10th the ripple, exactly.
        assert!(
            (large * 10.0 - small).abs() < 1e-9,
            "small={small} large={large}"
        );
    }

    #[test]
    fn higher_frequency_gives_less_ripple() {
        let mains_50 = capacitor_ripple_pp(0.5, 50.0, 1.0e-3).expect("valid");
        let mains_60 = capacitor_ripple_pp(0.5, 60.0, 1.0e-3).expect("valid");
        assert!(mains_60 < mains_50, "50hz={mains_50} 60hz={mains_60}");
    }

    #[test]
    fn full_wave_ripple_is_half_of_half_wave_for_same_mains() {
        // Same mains frequency: full-wave doubles the effective ripple
        // frequency, so its capacitor ripple is half the half-wave value.
        let mains = 60.0;
        let hw = capacitor_ripple_pp_for(Topology::HalfWave, 0.25, mains, 2.2e-3).expect("valid");
        let fw = capacitor_ripple_pp_for(Topology::FullWave, 0.25, mains, 2.2e-3).expect("valid");
        assert!((fw * 2.0 - hw).abs() < 1e-9, "hw={hw} fw={fw}");
    }

    #[test]
    fn ripple_multiplier_values() {
        assert!((Topology::HalfWave.ripple_frequency_multiplier() - 1.0).abs() < EPS);
        assert!((Topology::FullWave.ripple_frequency_multiplier() - 2.0).abs() < EPS);
    }

    #[test]
    fn zero_load_current_gives_zero_ripple() {
        let vr = capacitor_ripple_pp(0.0, 100.0, 1.0e-3).expect("valid");
        assert!(vr.abs() < EPS, "got {vr}");
    }

    // ---- Validation paths ------------------------------------------------

    #[test]
    fn rejects_non_positive_peak() {
        assert!(half_wave_vdc(0.0).is_err());
        assert!(full_wave_vdc(-1.0).is_err());
        assert!(half_wave_vrms(0.0).is_err());
        assert!(full_wave_vrms(-3.0).is_err());
    }

    #[test]
    fn rejects_non_finite_peak() {
        assert!(half_wave_vdc(f64::NAN).is_err());
        assert!(full_wave_vrms(f64::INFINITY).is_err());
    }

    #[test]
    fn capacitor_ripple_rejects_bad_inputs() {
        assert!(capacitor_ripple_pp(-1.0, 100.0, 1.0e-3).is_err()); // negative current
        assert!(capacitor_ripple_pp(1.0, 0.0, 1.0e-3).is_err()); // zero freq
        assert!(capacitor_ripple_pp(1.0, 100.0, 0.0).is_err()); // zero cap
        assert!(capacitor_ripple_pp(1.0, 100.0, f64::NAN).is_err()); // non-finite
    }

    #[test]
    fn ripple_factor_from_rejects_bad_inputs() {
        assert!(ripple_factor_from(1.0, 0.0).is_err()); // zero dc
        assert!(ripple_factor_from(-1.0, 1.0).is_err()); // negative rms
        assert!(ripple_factor_from(f64::INFINITY, 1.0).is_err()); // non-finite
    }

    #[test]
    fn topology_serde_roundtrip() {
        for t in [Topology::HalfWave, Topology::FullWave] {
            let json = serde_json::to_string(&t).expect("serialize");
            let back: Topology = serde_json::from_str(&json).expect("deserialize");
            assert_eq!(t, back);
        }
    }
}
