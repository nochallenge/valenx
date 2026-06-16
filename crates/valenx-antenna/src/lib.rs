//! # valenx-antenna — antenna-link fundamentals
//!
//! Closed-form models for the basic radio-link quantities that connect
//! an antenna's geometry and operating frequency to a received-power
//! estimate: the Friis transmission equation, aperture gain, decibel
//! conversions, the wavelength relation and a first-order beamwidth
//! estimate.
//!
//! ## What
//!
//! The [`wave`] module provides the `lambda = c / f` relation (and its
//! inverse), the medium-velocity form `lambda = v / f`, and the angular
//! wavenumber `k = 2*pi / lambda`; the exact speed of light
//! [`wave::SPEED_OF_LIGHT_M_S`] underpins the free-space forms.
//!
//! The [`gain`] module provides aperture gain `G = 4*pi*Ae / lambda^2`,
//! its inverse `Ae = G*lambda^2 / (4*pi)`, the effective aperture
//! `Ae = eta*A`, and the decibel-isotropic conversions
//! [`gain::to_dbi`] / [`gain::from_dbi`] (`G_dBi = 10*log10(G)`).
//!
//! The [`friis`] module provides the Friis transmission equation
//! `Pr/Pt = Gt*Gr*(lambda/(4*pi*d))^2`, received power, free-space
//! path loss in linear and dB form, and the inverse maximum-range solve
//! [`friis::max_range_m`] — the distance at which the received power
//! falls to a target receiver sensitivity.
//!
//! The [`beamwidth`] module provides the `theta ~ lambda / D` rule of
//! thumb, a `k`-weighted variant for common aperture tapers, and a
//! coarse directivity-from-beamwidth estimate.
//!
//! The [`spec`] module provides small serde-serializable
//! [`spec::LinkSpec`] / [`spec::LinkResult`] bundles that evaluate a
//! whole point-to-point link (including received power in dBm) in one
//! call.
//!
//! ## Model
//!
//! All relations are the standard far-field, free-space textbook closed
//! forms. With wavelength `lambda = c / f`:
//!
//! ```text
//! Friis:        Pr/Pt = Gt * Gr * (lambda / (4*pi*d))^2
//! aperture gain: G    = 4 * pi * Ae / lambda^2
//! decibels:     G_dBi = 10 * log10(G)
//! beamwidth:    theta ~ k * lambda / D
//! path loss:    FSPL  = (4*pi*d / lambda)^2
//! max range:    d_max = (lambda/(4*pi)) * sqrt(Pt*Gt*Gr / Pr_min)
//! ```
//!
//! Gains carried through the API are **linear** power ratios relative
//! to an isotropic radiator unless a function name says `dbi`. The
//! `1/d^2` distance roll-off and `1/f^2` (i.e. `lambda^2`) frequency
//! dependence of the Friis ratio, the `4*pi*Ae/lambda^2` aperture-gain
//! law, the `10*log10` decibel mapping and the "larger aperture means
//! higher gain and a narrower beam" trend are all exercised directly by
//! the test suite.
//!
//! ## Honest scope
//!
//! Research / educational grade. This crate implements the **standard
//! textbook closed-form / numerical** link-budget relations and nothing
//! more. It is **not** a clinical, medical, or production engineering
//! tool, and it is **not** a substitute for full-wave electromagnetic
//! simulation, measured antenna patterns, polarization / impedance
//! mismatch, atmospheric and multipath effects, near-field behaviour, or
//! regulatory RF-exposure and link analysis. The beamwidth and
//! directivity estimates in particular are deliberately first-order: a
//! real pattern depends on the full aperture illumination. Treat every
//! output as an order-of-magnitude teaching figure, not a certified
//! design value.
//!
//! ## Errors
//!
//! Every fallible function returns
//! [`Result<_, AntennaError>`](error::AntennaError); the error type
//! exposes stable [`code`](error::AntennaError::code) and
//! [`category`](error::AntennaError::category) accessors.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod beamwidth;
pub mod error;
pub mod friis;
pub mod gain;
pub mod spec;
pub mod wave;

pub use beamwidth::{
    beamwidth_k_rad, beamwidth_rad, directivity_from_beamwidths, radians_to_degrees,
    K_UNIFORM_CIRCULAR_NULL, K_UNIFORM_LINE_HPBW,
};
pub use error::{AntennaError, ErrorCategory};
pub use friis::{
    free_space_factor, free_space_path_loss, free_space_path_loss_db, max_range_m, power_ratio,
    received_power,
};
pub use gain::{aperture_from_gain, effective_aperture, from_dbi, gain_from_aperture, to_dbi};
pub use spec::{power_w_to_dbm, LinkResult, LinkSpec};
pub use wave::{
    frequency_from_wavelength, wavelength_from_frequency, wavelength_in_medium, wavenumber,
    SPEED_OF_LIGHT_M_S,
};

#[cfg(test)]
mod integration_tests {
    //! Cross-module ground-truth checks tying the relations together.

    use crate::beamwidth::beamwidth_rad;
    use crate::friis::{free_space_path_loss, power_ratio};
    use crate::gain::{aperture_from_gain, gain_from_aperture, to_dbi};
    use crate::wave::wavelength_from_frequency;

    #[test]
    fn aperture_gain_round_trips_through_dbi_and_back() {
        let lambda = wavelength_from_frequency(12.0e9).unwrap();
        let g = gain_from_aperture(2.0, lambda).unwrap();
        let g_dbi = to_dbi(g).unwrap();
        let g_lin = crate::gain::from_dbi(g_dbi).unwrap();
        let ae = aperture_from_gain(g_lin, lambda).unwrap();
        assert!((ae - 2.0).abs() < 1e-9, "aperture round-trip failed: {ae}");
    }

    #[test]
    fn friis_equals_combined_gain_over_path_loss() {
        // Pr/Pt = Gt*Gr / FSPL_linear — a restatement of Friis that the
        // two independent code paths must agree on.
        let lambda = wavelength_from_frequency(5.8e9).unwrap();
        let d = 2_500.0;
        let (gt, gr) = (6.3, 12.0);
        let ratio = power_ratio(gt, gr, lambda, d).unwrap();
        let fspl = free_space_path_loss(lambda, d).unwrap();
        let via_loss = gt * gr / fspl;
        assert!(
            (ratio - via_loss).abs() / ratio < 1e-12,
            "Friis vs path-loss mismatch: {ratio} vs {via_loss}"
        );
    }

    #[test]
    fn larger_aperture_higher_gain_and_narrower_beam_together() {
        // The headline qualitative claim: a bigger aperture both raises
        // gain and narrows the beam at fixed wavelength.
        let lambda = 0.05;
        let small_d = 0.5;
        let large_d = 5.0;

        // Use a square aperture so Ae scales as D^2 for the gain side.
        let g_small = gain_from_aperture(small_d * small_d, lambda).unwrap();
        let g_large = gain_from_aperture(large_d * large_d, lambda).unwrap();
        assert!(g_large > g_small, "gain should increase with aperture");

        let beam_small = beamwidth_rad(lambda, small_d).unwrap();
        let beam_large = beamwidth_rad(lambda, large_d).unwrap();
        assert!(
            beam_large < beam_small,
            "beam should narrow with larger aperture"
        );
    }
}
