//! The transformer EMF equation — RMS induced voltage from core flux.
//!
//! When a winding of `N` turns links a sinusoidal mutual flux of peak
//! value `Phi_max` (webers) oscillating at frequency `f` (hertz), the
//! RMS induced electromotive force is
//!
//! ```text
//! E_rms = (2*pi / sqrt(2)) * f * N * Phi_max = sqrt(2) * pi * f * N * Phi_max
//!       ~= 4.44 * f * N * Phi_max
//! ```
//!
//! the classic **transformer EMF equation**. The leading constant is
//! exactly `sqrt(2)*pi = 4.442_882_9...` ([`FORM_FACTOR`]); the familiar
//! `4.44` of the textbooks is just that value rounded to three figures,
//! so this module evaluates the *exact* form. It arises as `(2*pi)/sqrt(2)`:
//! the `2*pi` turns frequency into the angular frequency of the peak EMF
//! `N*Phi_max*omega`, and the `1/sqrt(2)` turns that peak into an RMS value.
//!
//! ## Why it grounds the turns ratio
//!
//! Both windings of a transformer link the *same* mutual flux at the
//! *same* frequency, so the EMF *per turn* — [`volts_per_turn`],
//! `E/N = sqrt(2)*pi*f*Phi_max` — is identical on the primary and the
//! secondary. Dividing the two winding EMFs cancels that common factor:
//!
//! ```text
//! Ep / Es = (sqrt(2) pi f Np Phi) / (sqrt(2) pi f Ns Phi) = Np / Ns = a
//! ```
//!
//! which is exactly the ideal voltage law of [`crate::ratio`]. The
//! axiomatic turns-ratio relation there is the consequence of this
//! equation.
//!
//! ## Honest scope
//!
//! This is the ideal sinusoidal-excitation EMF law: a single mutual
//! flux, one frequency, a perfectly sinusoidal flux waveform, and no
//! saturation, leakage, or core loss. `Phi_max` is the *peak* flux
//! through the core, not a flux density; multiply a peak flux density
//! `B_max` (teslas) by the core cross-sectional area (square metres) to
//! obtain it. As with the rest of the crate this is a teaching / sizing
//! aid, not a core-design tool.

use crate::error::TransformerError;

/// The exact EMF-equation form factor `sqrt(2)*pi = 4.442_882_938...`.
///
/// This is the leading constant of the transformer EMF equation
/// `E_rms = FORM_FACTOR * f * N * Phi_max`. It is `(2*pi)/sqrt(2)`: the
/// `2*pi` converts frequency to the angular frequency of the peak EMF
/// `N*Phi_max*omega`, and the `1/sqrt(2)` converts that peak to an RMS
/// value. The textbook constant `4.44` is this number rounded to three
/// significant figures.
pub const FORM_FACTOR: f64 = std::f64::consts::SQRT_2 * std::f64::consts::PI;

/// RMS induced EMF of a winding from the transformer EMF equation,
/// `E_rms = sqrt(2)*pi*f*N*Phi_max`, in volts.
///
/// `frequency_hz` is the excitation frequency `f` (Hz), `turns` the
/// winding turn count `N`, and `peak_flux_wb` the peak mutual flux
/// `Phi_max` through the core (webers). Because a weber is a
/// volt-second, `Hz * Wb = V`, so the result is in volts.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if `frequency_hz` or `turns` is
/// not finite and strictly positive, or if `peak_flux_wb` is not finite
/// and non-negative.
pub fn induced_emf_rms(
    frequency_hz: f64,
    turns: f64,
    peak_flux_wb: f64,
) -> Result<f64, TransformerError> {
    if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
        return Err(TransformerError::invalid(
            "frequency_hz",
            format!("frequency must be finite and positive, got {frequency_hz}"),
        ));
    }
    if !turns.is_finite() || turns <= 0.0 {
        return Err(TransformerError::invalid(
            "turns",
            format!("number of turns must be finite and positive, got {turns}"),
        ));
    }
    if !peak_flux_wb.is_finite() || peak_flux_wb < 0.0 {
        return Err(TransformerError::invalid(
            "peak_flux_wb",
            format!("peak flux must be finite and non-negative, got {peak_flux_wb}"),
        ));
    }
    Ok(FORM_FACTOR * frequency_hz * turns * peak_flux_wb)
}

/// The number of turns needed to induce a target RMS EMF at a given
/// frequency and peak core flux, inverting [`induced_emf_rms`]:
/// `N = E_rms / (sqrt(2)*pi*f*Phi_max)`.
///
/// This is the design form of the EMF equation: pick the operating
/// frequency and the peak core flux (set by the core area and the
/// saturation limit of the material), and read off how many turns a
/// winding needs to develop the wanted RMS voltage. The result is a real
/// turn count and is generally not an integer; a real winding rounds it.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if `emf_rms` is not finite and
/// non-negative, if `frequency_hz` is not finite and strictly positive,
/// or if `peak_flux_wb` is not finite and strictly positive (it is the
/// ratio denominator).
pub fn turns_for_emf(
    emf_rms: f64,
    frequency_hz: f64,
    peak_flux_wb: f64,
) -> Result<f64, TransformerError> {
    if !emf_rms.is_finite() || emf_rms < 0.0 {
        return Err(TransformerError::invalid(
            "emf_rms",
            format!("EMF must be finite and non-negative, got {emf_rms}"),
        ));
    }
    if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
        return Err(TransformerError::invalid(
            "frequency_hz",
            format!("frequency must be finite and positive, got {frequency_hz}"),
        ));
    }
    if !peak_flux_wb.is_finite() || peak_flux_wb <= 0.0 {
        return Err(TransformerError::invalid(
            "peak_flux_wb",
            format!("peak flux must be finite and positive, got {peak_flux_wb}"),
        ));
    }
    Ok(emf_rms / (FORM_FACTOR * frequency_hz * peak_flux_wb))
}

/// The peak mutual flux a winding drives to develop a given RMS EMF at a
/// given frequency, inverting [`induced_emf_rms`]:
/// `Phi_max = E_rms / (sqrt(2)*pi*f*N)`, in webers.
///
/// Use this to check a design against core saturation: divide the
/// returned peak flux by the core cross-sectional area to get the peak
/// flux density `B_max`, and compare it with the saturation flux density
/// of the core material.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if `emf_rms` is not finite and
/// non-negative, if `frequency_hz` is not finite and strictly positive,
/// or if `turns` is not finite and strictly positive (it is the ratio
/// denominator).
pub fn peak_flux_for_emf(
    emf_rms: f64,
    frequency_hz: f64,
    turns: f64,
) -> Result<f64, TransformerError> {
    if !emf_rms.is_finite() || emf_rms < 0.0 {
        return Err(TransformerError::invalid(
            "emf_rms",
            format!("EMF must be finite and non-negative, got {emf_rms}"),
        ));
    }
    if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
        return Err(TransformerError::invalid(
            "frequency_hz",
            format!("frequency must be finite and positive, got {frequency_hz}"),
        ));
    }
    if !turns.is_finite() || turns <= 0.0 {
        return Err(TransformerError::invalid(
            "turns",
            format!("number of turns must be finite and positive, got {turns}"),
        ));
    }
    Ok(emf_rms / (FORM_FACTOR * frequency_hz * turns))
}

/// The induced EMF *per turn*, `E/N = sqrt(2)*pi*f*Phi_max`, in volts per
/// turn.
///
/// Because every winding on the same core links the same mutual flux
/// `Phi_max` at the same frequency `f`, this volts-per-turn figure is
/// identical for the primary and the secondary — which is precisely why
/// `Vp/Vs = Np/Ns = a` in [`crate::ratio`]. It is also the natural
/// design knob: multiply by a winding's turn count to get that winding's
/// RMS voltage.
///
/// # Errors
///
/// Returns [`TransformerError::Invalid`] if `frequency_hz` is not finite
/// and strictly positive, or if `peak_flux_wb` is not finite and
/// non-negative.
pub fn volts_per_turn(frequency_hz: f64, peak_flux_wb: f64) -> Result<f64, TransformerError> {
    if !frequency_hz.is_finite() || frequency_hz <= 0.0 {
        return Err(TransformerError::invalid(
            "frequency_hz",
            format!("frequency must be finite and positive, got {frequency_hz}"),
        ));
    }
    if !peak_flux_wb.is_finite() || peak_flux_wb < 0.0 {
        return Err(TransformerError::invalid(
            "peak_flux_wb",
            format!("peak flux must be finite and non-negative, got {peak_flux_wb}"),
        ));
    }
    Ok(FORM_FACTOR * frequency_hz * peak_flux_wb)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ratio::TurnsRatio;

    /// Relative tolerance for the analytic float checks.
    const EPS: f64 = 1e-9;

    #[test]
    fn form_factor_is_root_two_pi_near_textbook_four_point_four_four() {
        // The exact constant is sqrt(2)*pi; the textbook 4.44 rounds it.
        // Pinning it to 4.4428829... within 1e-12 also confirms it lies
        // just above the rounded 4.44 of the textbooks.
        assert!((FORM_FACTOR - 4.442_882_938_158_366).abs() < 1e-12);
        // Cross-check the constant against an independent (2*pi)/sqrt(2)
        // construction (a runtime product, not the const, so this is a
        // real comparison rather than a constant assertion).
        let two_pi_over_root_two = (2.0 * std::f64::consts::PI) / 2.0_f64.sqrt();
        assert!((FORM_FACTOR - two_pi_over_root_two).abs() < 1e-12);
    }

    #[test]
    fn emf_equation_closed_form() {
        // f = 50 Hz, N = 100 turns, Phi_max = 0.01 Wb.
        // E = sqrt(2)*pi*50*100*0.01 = 222.1441... V; the rounded
        // textbook constant 4.44 gives 222.0 V, so this is just above it.
        let e = induced_emf_rms(50.0, 100.0, 0.01).unwrap();
        assert!((e - 222.144_146_907_918_3).abs() < 1e-6, "E got {e}");
        assert!(
            (e - 4.44 * 50.0 * 100.0 * 0.01).abs() < 0.2,
            "near 4.44 form"
        );
    }

    #[test]
    fn emf_is_linear_in_each_factor() {
        let base = induced_emf_rms(50.0, 100.0, 0.01).unwrap();
        // Doubling the frequency doubles the EMF.
        let f2 = induced_emf_rms(100.0, 100.0, 0.01).unwrap();
        assert!(
            (f2 / base - 2.0).abs() < EPS,
            "freq scaling got {}",
            f2 / base
        );
        // Doubling the turns doubles the EMF.
        let n2 = induced_emf_rms(50.0, 200.0, 0.01).unwrap();
        assert!(
            (n2 / base - 2.0).abs() < EPS,
            "turns scaling got {}",
            n2 / base
        );
        // Doubling the flux doubles the EMF.
        let p2 = induced_emf_rms(50.0, 100.0, 0.02).unwrap();
        assert!(
            (p2 / base - 2.0).abs() < EPS,
            "flux scaling got {}",
            p2 / base
        );
    }

    #[test]
    fn turns_for_emf_inverts_induced_emf() {
        // Round-trip: turns -> EMF -> turns recovers the turn count.
        let (f, n, flux) = (60.0, 275.0, 0.008);
        let e = induced_emf_rms(f, n, flux).unwrap();
        let n_back = turns_for_emf(e, f, flux).unwrap();
        assert!((n_back - n).abs() < 1e-9, "turns round-trip got {n_back}");
    }

    #[test]
    fn peak_flux_for_emf_inverts_induced_emf() {
        // Round-trip: flux -> EMF -> flux recovers the peak flux.
        let (f, n, flux) = (400.0, 32.0, 0.0015);
        let e = induced_emf_rms(f, n, flux).unwrap();
        let flux_back = peak_flux_for_emf(e, f, n).unwrap();
        assert!(
            (flux_back - flux).abs() < 1e-12,
            "flux round-trip got {flux_back}"
        );
    }

    #[test]
    fn volts_per_turn_times_turns_is_the_emf() {
        let (f, n, flux) = (50.0, 100.0, 0.01);
        let vpt = volts_per_turn(f, flux).unwrap();
        let e = induced_emf_rms(f, n, flux).unwrap();
        assert!((vpt * n - e).abs() < 1e-9, "E/N * N != E: {}", vpt * n);
    }

    #[test]
    fn emf_equation_reproduces_the_turns_ratio_law() {
        // GOLD derivable identity: both windings link the same flux at
        // the same frequency, so the per-turn EMF is shared and the EMF
        // ratio equals the turns ratio of crate::ratio.
        let (np, ns) = (240.0, 24.0); // a = 10, step-down
        let (f, flux) = (60.0, 0.005);
        let ep = induced_emf_rms(f, np, flux).unwrap();
        let es = induced_emf_rms(f, ns, flux).unwrap();

        let a = TurnsRatio::from_turns(np, ns).unwrap();
        // Ep / Es == Np / Ns == a.
        assert!((ep / es - a.ratio()).abs() < 1e-9, "Ep/Es got {}", ep / es);
        // The ratio module's voltage law maps Ep straight onto Es.
        let es_via_ratio = a.secondary_voltage(ep).unwrap();
        assert!(
            (es_via_ratio - es).abs() < 1e-9,
            "Es via ratio got {es_via_ratio}"
        );

        // The volts-per-turn is identical on both windings.
        let vpt = volts_per_turn(f, flux).unwrap();
        assert!((ep / np - vpt).abs() < 1e-9, "primary V/turn");
        assert!((es / ns - vpt).abs() < 1e-9, "secondary V/turn");
    }

    #[test]
    fn zero_flux_gives_zero_emf() {
        // A non-negative-magnitude convention: zero peak flux is allowed
        // on the forward path and yields zero EMF.
        assert!(induced_emf_rms(50.0, 100.0, 0.0).unwrap().abs() < 1e-12);
        assert!(volts_per_turn(50.0, 0.0).unwrap().abs() < 1e-12);
    }

    #[test]
    fn rejects_out_of_domain_inputs() {
        // Forward law: f and N strictly positive, flux non-negative.
        assert!(induced_emf_rms(0.0, 100.0, 0.01).is_err());
        assert!(induced_emf_rms(50.0, 0.0, 0.01).is_err());
        assert!(induced_emf_rms(50.0, 100.0, -0.01).is_err());
        assert!(induced_emf_rms(f64::NAN, 100.0, 0.01).is_err());
        assert!(induced_emf_rms(f64::INFINITY, 100.0, 0.01).is_err());

        // Inverses reject a zero denominator (flux / turns).
        assert!(turns_for_emf(100.0, 50.0, 0.0).is_err());
        assert!(turns_for_emf(-1.0, 50.0, 0.01).is_err());
        assert!(peak_flux_for_emf(100.0, 50.0, 0.0).is_err());
        assert!(peak_flux_for_emf(100.0, 0.0, 10.0).is_err());

        assert!(volts_per_turn(-1.0, 0.01).is_err());
        assert!(volts_per_turn(50.0, f64::NAN).is_err());
    }
}
