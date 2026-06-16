//! Gain-bandwidth relations for a single-pole op-amp.
//!
//! A real op-amp's open-loop response rolls off as a single dominant
//! pole, so the product of closed-loop gain and closed-loop bandwidth is
//! a constant: the **gain-bandwidth product** (GBW, sometimes the
//! unity-gain frequency `f_T`). Three textbook identities follow:
//!
//! 1. `GBW = gain · bandwidth`
//! 2. closed-loop bandwidth `= GBW / |gain|`
//! 3. unity-gain bandwidth `= GBW` (the gain-1 case)
//!
//! These are first-order approximations. Real parts add extra poles,
//! slew-rate limiting and phase margin effects this model ignores.

use crate::error::{ensure_positive, Result};

/// The gain-bandwidth product of a single-pole op-amp.
///
/// Stored in hertz (the unity-gain frequency). Construct with
/// [`Gbw::new`] or derive one from a known gain / bandwidth pair via
/// [`Gbw::from_gain_bandwidth`].
#[derive(Copy, Clone, Debug, PartialEq, PartialOrd)]
pub struct Gbw {
    hz: f64,
}

impl Gbw {
    /// Construct from a gain-bandwidth product in hertz.
    ///
    /// # Errors
    ///
    /// Returns [`OpAmpError`](crate::OpAmpError) if `hz` is non-finite
    /// or `<= 0`.
    pub fn new(hz: f64) -> Result<Self> {
        Ok(Self {
            hz: ensure_positive("gbw_hz", hz)?,
        })
    }

    /// Derive the gain-bandwidth product from a measured closed-loop
    /// gain magnitude and its bandwidth: `GBW = |gain| · bandwidth`.
    ///
    /// # Errors
    ///
    /// Returns [`OpAmpError`](crate::OpAmpError) if either argument is
    /// non-finite or `<= 0`.
    pub fn from_gain_bandwidth(gain_magnitude: f64, bandwidth_hz: f64) -> Result<Self> {
        let g = ensure_positive("gain_magnitude", gain_magnitude)?;
        let bw = ensure_positive("bandwidth_hz", bandwidth_hz)?;
        Self::new(g * bw)
    }

    /// The gain-bandwidth product in hertz.
    pub fn hz(&self) -> f64 {
        self.hz
    }

    /// Closed-loop bandwidth for a given gain magnitude:
    /// `BW = GBW / |gain|`.
    ///
    /// Higher closed-loop gain trades away bandwidth in exact inverse
    /// proportion.
    ///
    /// # Errors
    ///
    /// Returns [`OpAmpError`](crate::OpAmpError) if `gain_magnitude` is
    /// non-finite or `<= 0`.
    pub fn closed_loop_bandwidth(&self, gain_magnitude: f64) -> Result<f64> {
        let g = ensure_positive("gain_magnitude", gain_magnitude)?;
        Ok(self.hz / g)
    }

    /// The maximum closed-loop gain magnitude that still achieves a
    /// required closed-loop bandwidth: `|gain| = GBW / bandwidth`.
    ///
    /// This is the exact inverse of
    /// [`closed_loop_bandwidth`](Self::closed_loop_bandwidth): where that
    /// trades gain away for bandwidth, this answers the dual design
    /// question — "given a bandwidth I must hit, how much gain can I
    /// afford?". Demanding more bandwidth lowers the affordable gain in
    /// exact inverse proportion; at `bandwidth = GBW` (the unity-gain
    /// frequency) the gain falls to `1`, and a bandwidth beyond `GBW`
    /// implies an attenuating (`gain < 1`) configuration, mirroring how
    /// [`closed_loop_bandwidth`](Self::closed_loop_bandwidth) accepts a
    /// sub-unity gain.
    ///
    /// # Errors
    ///
    /// Returns [`OpAmpError`](crate::OpAmpError) if `bandwidth_hz` is
    /// non-finite or `<= 0`.
    pub fn max_gain_for_bandwidth(&self, bandwidth_hz: f64) -> Result<f64> {
        let bw = ensure_positive("bandwidth_hz", bandwidth_hz)?;
        Ok(self.hz / bw)
    }

    /// Unity-gain bandwidth — the frequency at which the open-loop gain
    /// falls to `1`. For a single-pole op-amp this equals the
    /// gain-bandwidth product itself.
    pub fn unity_gain_bandwidth(&self) -> f64 {
        self.hz
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    #[test]
    fn max_gain_hand_value() {
        // GBW = 1 MHz, need 10 kHz bandwidth -> max gain 100.
        let gbw = Gbw::new(1.0e6).unwrap();
        let g = gbw.max_gain_for_bandwidth(10_000.0).unwrap();
        assert!((g - 100.0).abs() < EPS, "gain = {g}");
    }

    #[test]
    fn max_gain_inverts_closed_loop_bandwidth_both_ways() {
        let gbw = Gbw::new(2.5e6).unwrap();
        // gain -> bandwidth -> gain.
        let g0 = 47.0;
        let bw = gbw.closed_loop_bandwidth(g0).unwrap();
        let g_back = gbw.max_gain_for_bandwidth(bw).unwrap();
        assert!((g_back - g0).abs() < 1e-9 * g0, "gain round-trip {g_back}");
        // bandwidth -> gain -> bandwidth.
        let bw0 = 31_250.0;
        let g = gbw.max_gain_for_bandwidth(bw0).unwrap();
        let bw_back = gbw.closed_loop_bandwidth(g).unwrap();
        assert!(
            (bw_back - bw0).abs() < 1e-6 * bw0,
            "bandwidth round-trip {bw_back}"
        );
    }

    #[test]
    fn gain_times_bandwidth_is_gbw() {
        // The defining identity: |gain| * bandwidth == GBW.
        let gbw = Gbw::new(4.0e6).unwrap();
        for &bw in &[1.0e3_f64, 1.0e4, 2.0e5, 1.0e6] {
            let g = gbw.max_gain_for_bandwidth(bw).unwrap();
            assert!((g * bw - gbw.hz()).abs() < 1.0, "gain*bw != GBW at bw={bw}");
        }
    }

    #[test]
    fn bandwidth_equal_to_gbw_gives_unity_gain() {
        let gbw = Gbw::new(1.5e6).unwrap();
        let g = gbw
            .max_gain_for_bandwidth(gbw.unity_gain_bandwidth())
            .unwrap();
        assert!(
            (g - 1.0).abs() < EPS,
            "unity-gain bandwidth should give gain 1, got {g}"
        );
    }

    #[test]
    fn more_bandwidth_means_less_gain() {
        let gbw = Gbw::new(1.0e6).unwrap();
        let lo_bw = gbw.max_gain_for_bandwidth(1.0e4).unwrap();
        let hi_bw = gbw.max_gain_for_bandwidth(1.0e5).unwrap();
        assert!(
            hi_bw < lo_bw,
            "more bandwidth should allow less gain: {hi_bw} vs {lo_bw}"
        );
    }

    #[test]
    fn max_gain_rejects_bad_bandwidth() {
        let gbw = Gbw::new(1.0e6).unwrap();
        assert!(gbw.max_gain_for_bandwidth(0.0).is_err());
        assert!(gbw.max_gain_for_bandwidth(-1.0).is_err());
        assert!(gbw.max_gain_for_bandwidth(f64::NAN).is_err());
        assert!(gbw.max_gain_for_bandwidth(f64::INFINITY).is_err());
    }
}
