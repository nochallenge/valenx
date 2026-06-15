//! Serializable parameter bundles for a point-to-point antenna link.
//!
//! These small structs let a host (a GUI panel, an LLM/MCP request
//! surface, a saved scenario file) round-trip a link description and
//! evaluate it in one call via [`LinkSpec::evaluate`]. They are thin
//! wrappers over the free functions in [`crate::friis`], [`crate::gain`]
//! and [`crate::wave`].

use crate::error::{require_positive, AntennaError};
use crate::friis::{free_space_path_loss_db, power_ratio, received_power};
use crate::gain::to_dbi;
use crate::wave::wavelength_from_frequency;
use serde::{Deserialize, Serialize};

/// A complete free-space point-to-point link description.
///
/// All gains are given as **linear** ratios; convert from dBi with
/// [`crate::gain::from_dbi`] before constructing.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct LinkSpec {
    /// Operating frequency in hertz.
    pub freq_hz: f64,
    /// Transmit power in watts.
    pub tx_power_w: f64,
    /// Transmit antenna gain (linear ratio).
    pub tx_gain: f64,
    /// Receive antenna gain (linear ratio).
    pub rx_gain: f64,
    /// Line-of-sight separation in metres.
    pub distance_m: f64,
}

/// The evaluated result of a [`LinkSpec`].
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct LinkResult {
    /// Free-space wavelength in metres (`lambda = c / f`).
    pub wavelength_m: f64,
    /// Friis power ratio `Pr/Pt` (dimensionless).
    pub power_ratio: f64,
    /// Received power in watts.
    pub rx_power_w: f64,
    /// Received power in dBm (decibels relative to one milliwatt).
    pub rx_power_dbm: f64,
    /// Free-space path loss in dB.
    pub path_loss_db: f64,
}

impl LinkSpec {
    /// Evaluate the link, returning the derived wavelength, Friis power
    /// ratio, received power (watts and dBm) and free-space path loss.
    ///
    /// # Errors
    ///
    /// Returns an error if any field is out of domain (non-positive
    /// frequency / power / distance, or non-finite / negative gain),
    /// propagated from the underlying closed-form functions.
    pub fn evaluate(&self) -> Result<LinkResult, AntennaError> {
        let wavelength_m = wavelength_from_frequency(self.freq_hz)?;
        let power_ratio = power_ratio(self.tx_gain, self.rx_gain, wavelength_m, self.distance_m)?;
        let rx_power_w = received_power(
            self.tx_power_w,
            self.tx_gain,
            self.rx_gain,
            wavelength_m,
            self.distance_m,
        )?;
        // rx_power_w is strictly positive here (tx_power_w > 0 and the
        // ratio is > 0 because the gains were validated > 0 implicitly
        // through the path; guard anyway for non-positive ratio with a
        // zero gain).
        let rx_power_dbm = power_w_to_dbm(rx_power_w)?;
        let path_loss_db = free_space_path_loss_db(wavelength_m, self.distance_m)?;
        Ok(LinkResult {
            wavelength_m,
            power_ratio,
            rx_power_w,
            rx_power_dbm,
            path_loss_db,
        })
    }
}

/// Convert an absolute power in watts to dBm:
/// `P_dBm = 10*log10(P_W / 1e-3) = 10*log10(P_W) + 30`.
///
/// # Errors
///
/// Returns an error if `power_w` is not finite and strictly positive.
pub fn power_w_to_dbm(power_w: f64) -> Result<f64, AntennaError> {
    let p = require_positive("power_w", power_w)?;
    // The gain helper computes 10*log10 with the same positivity gate;
    // reuse it for the dB form and add the +30 dBm offset.
    Ok(to_dbi(p)? + 30.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn watts_to_dbm_known_values() {
        // 1 W = 30 dBm, 1 mW = 0 dBm, 1 microwatt = -30 dBm.
        assert!((power_w_to_dbm(1.0).unwrap() - 30.0).abs() < 1e-9);
        assert!((power_w_to_dbm(1.0e-3).unwrap() - 0.0).abs() < 1e-9);
        assert!((power_w_to_dbm(1.0e-6).unwrap() - (-30.0)).abs() < 1e-9);
    }

    #[test]
    fn evaluate_matches_components() {
        let spec = LinkSpec {
            freq_hz: 2.4e9,
            tx_power_w: 0.1,
            tx_gain: 2.0,
            rx_gain: 2.0,
            distance_m: 100.0,
        };
        let r = spec.evaluate().unwrap();

        let lambda = wavelength_from_frequency(2.4e9).unwrap();
        assert!((r.wavelength_m - lambda).abs() < 1e-12);

        let expect_ratio = power_ratio(2.0, 2.0, lambda, 100.0).unwrap();
        assert!((r.power_ratio - expect_ratio).abs() < 1e-15);

        let expect_pr = 0.1 * expect_ratio;
        assert!((r.rx_power_w - expect_pr).abs() < 1e-18);

        // Path loss should be ~80.05 dB (2.4 GHz, 100 m).
        assert!(
            (r.path_loss_db - 80.05).abs() < 0.05,
            "got {}",
            r.path_loss_db
        );
    }

    #[test]
    fn evaluate_rejects_bad_spec() {
        let bad = LinkSpec {
            freq_hz: -1.0,
            tx_power_w: 1.0,
            tx_gain: 1.0,
            rx_gain: 1.0,
            distance_m: 10.0,
        };
        assert!(bad.evaluate().is_err());
    }

    #[test]
    fn spec_serde_roundtrip() {
        let spec = LinkSpec {
            freq_hz: 1.0e9,
            tx_power_w: 5.0,
            tx_gain: 10.0,
            rx_gain: 3.0,
            distance_m: 1_000.0,
        };
        let json = serde_json::to_string(&spec).unwrap();
        let back: LinkSpec = serde_json::from_str(&json).unwrap();
        assert_eq!(spec, back);
    }
}
