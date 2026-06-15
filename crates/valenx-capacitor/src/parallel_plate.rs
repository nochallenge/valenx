//! Parallel-plate capacitance and stored electrostatic energy.
//!
//! ## Model
//!
//! For two conductive plates of overlapping area `A` (m^2) separated by a
//! gap `d` (m) filled with a linear dielectric of relative permittivity
//! `eps_r`, the idealised capacitance is
//!
//! ```text
//! C = eps0 * eps_r * A / d
//! ```
//!
//! where `eps0` is the [vacuum permittivity](VACUUM_PERMITTIVITY). The
//! energy stored at terminal voltage `V` (volts) is
//!
//! ```text
//! E = 1/2 * C * V^2
//! ```
//!
//! ## Honest scope
//!
//! This is the textbook infinite-parallel-plate result. It assumes a
//! uniform field confined between the plates and therefore ignores
//! fringing at the plate edges (which makes a real capacitor slightly
//! larger than this formula predicts), dielectric non-linearity, loss and
//! breakdown. It is suitable for teaching and order-of-magnitude work, not
//! for production capacitor design.

use crate::error::{CapacitorError, Result};

/// Vacuum permittivity (electric constant) `eps0`, in farads per metre
/// (F/m). 2018 CODATA value.
pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_812_8e-12;

/// Capacitance of an ideal parallel-plate capacitor, in farads (F).
///
/// Computes `C = eps0 * eps_r * A / d`.
///
/// # Parameters
///
/// - `eps_r` — relative permittivity of the dielectric (dimensionless,
///   `1.0` for vacuum / dry air). Must be `>= 1`.
/// - `area_m2` — overlapping plate area `A`, in m^2. Must be `> 0`.
/// - `gap_m` — plate separation `d`, in metres. Must be `> 0`.
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if `eps_r < 1`, or if
/// `area_m2` or `gap_m` is not strictly positive and finite.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::parallel_plate::capacitance;
///
/// // 1 cm^2 plates, 1 mm apart, in vacuum.
/// let c = capacitance(1.0, 1.0e-4, 1.0e-3).unwrap();
/// assert!((c - 8.854_187_812_8e-13).abs() < 1e-24);
/// ```
pub fn capacitance(eps_r: f64, area_m2: f64, gap_m: f64) -> Result<f64> {
    if !(eps_r.is_finite() && eps_r >= 1.0) {
        return Err(CapacitorError::InvalidParameter {
            name: "eps_r",
            value: eps_r,
            reason: "relative permittivity must be a finite value >= 1",
        });
    }
    let area = CapacitorError::require_positive("area_m2", area_m2, "plate area must be positive")?;
    let gap = CapacitorError::require_positive("gap_m", gap_m, "plate gap must be positive")?;
    Ok(VACUUM_PERMITTIVITY * eps_r * area / gap)
}

/// Electrostatic energy stored in a capacitor, in joules (J).
///
/// Computes `E = 1/2 * C * V^2`.
///
/// # Parameters
///
/// - `capacitance_f` — capacitance `C`, in farads. Must be `> 0`.
/// - `voltage_v` — terminal voltage `V`, in volts. Any finite value is
///   accepted; energy depends on `V^2`, so polarity does not matter.
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if `capacitance_f` is not
/// strictly positive and finite, or if `voltage_v` is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::parallel_plate::stored_energy;
///
/// // 100 uF charged to 10 V stores 5 mJ.
/// let e = stored_energy(100.0e-6, 10.0).unwrap();
/// assert!((e - 5.0e-3).abs() < 1e-12);
/// ```
pub fn stored_energy(capacitance_f: f64, voltage_v: f64) -> Result<f64> {
    let c = CapacitorError::require_positive(
        "capacitance_f",
        capacitance_f,
        "capacitance must be positive",
    )?;
    if !voltage_v.is_finite() {
        return Err(CapacitorError::InvalidParameter {
            name: "voltage_v",
            value: voltage_v,
            reason: "voltage must be finite",
        });
    }
    Ok(0.5 * c * voltage_v * voltage_v)
}

/// Charge stored on a capacitor at a given voltage, in coulombs (C).
///
/// Computes the defining relation `Q = C * V`.
///
/// # Parameters
///
/// - `capacitance_f` — capacitance `C`, in farads. Must be `> 0`.
/// - `voltage_v` — terminal voltage `V`, in volts. Any finite value is
///   accepted.
///
/// # Errors
///
/// Returns [`CapacitorError::InvalidParameter`] if `capacitance_f` is not
/// strictly positive and finite, or if `voltage_v` is non-finite.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::parallel_plate::charge;
///
/// // 1 uF at 5 V holds 5 uC.
/// let q = charge(1.0e-6, 5.0).unwrap();
/// assert!((q - 5.0e-6).abs() < 1e-18);
/// ```
pub fn charge(capacitance_f: f64, voltage_v: f64) -> Result<f64> {
    let c = CapacitorError::require_positive(
        "capacitance_f",
        capacitance_f,
        "capacitance must be positive",
    )?;
    if !voltage_v.is_finite() {
        return Err(CapacitorError::InvalidParameter {
            name: "voltage_v",
            value: voltage_v,
            reason: "voltage must be finite",
        });
    }
    Ok(c * voltage_v)
}
