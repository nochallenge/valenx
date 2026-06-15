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

    /// Unity-gain bandwidth — the frequency at which the open-loop gain
    /// falls to `1`. For a single-pole op-amp this equals the
    /// gain-bandwidth product itself.
    pub fn unity_gain_bandwidth(&self) -> f64 {
        self.hz
    }
}
