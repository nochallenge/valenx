//! Series and parallel combination of ideal capacitors.
//!
//! ## Model
//!
//! Capacitors combine oppositely to resistors. For `n` capacitors:
//!
//! ```text
//! parallel:  C_total = C_1 + C_2 + ... + C_n
//! series:    1 / C_total = 1/C_1 + 1/C_2 + ... + 1/C_n
//! ```
//!
//! Parallel capacitances add (the plate areas effectively sum); series
//! capacitances combine reciprocally, so the total is always *smaller*
//! than the smallest branch.
//!
//! ## Honest scope
//!
//! These are the ideal network identities. They assume every capacitor is
//! a pure capacitance with no tolerance spread, no leakage and no series
//! parasitics, and they say nothing about the voltage rating of the
//! combination (in a real series stack the applied voltage divides
//! unevenly across mismatched parts).

use crate::error::{CapacitorError, Result};

/// Equivalent capacitance of capacitors wired in **parallel**, in farads.
///
/// Computes `C_total = C_1 + C_2 + ... + C_n`.
///
/// # Parameters
///
/// - `caps_f` — the branch capacitances, in farads. The slice must be
///   non-empty and every entry must be `> 0`.
///
/// # Errors
///
/// Returns [`CapacitorError::EmptyNetwork`] if `caps_f` is empty, or
/// [`CapacitorError::InvalidParameter`] if any entry is not strictly
/// positive and finite.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::network::parallel;
///
/// // 1 uF || 2 uF || 3 uF = 6 uF.
/// let c = parallel(&[1.0e-6, 2.0e-6, 3.0e-6]).unwrap();
/// assert!((c - 6.0e-6).abs() < 1e-18);
/// ```
pub fn parallel(caps_f: &[f64]) -> Result<f64> {
    if caps_f.is_empty() {
        return Err(CapacitorError::EmptyNetwork(
            "parallel combination needs at least one capacitor",
        ));
    }
    let mut total = 0.0;
    for (i, &c) in caps_f.iter().enumerate() {
        if !(c.is_finite() && c > 0.0) {
            return Err(CapacitorError::InvalidParameter {
                name: index_name(i),
                value: c,
                reason: "each branch capacitance must be positive",
            });
        }
        total += c;
    }
    Ok(total)
}

/// Equivalent capacitance of capacitors wired in **series**, in farads.
///
/// Computes `1 / C_total = 1/C_1 + ... + 1/C_n`, i.e. the reciprocal of
/// the sum of reciprocals. The result is always smaller than the smallest
/// branch.
///
/// # Parameters
///
/// - `caps_f` — the branch capacitances, in farads. The slice must be
///   non-empty and every entry must be `> 0`.
///
/// # Errors
///
/// Returns [`CapacitorError::EmptyNetwork`] if `caps_f` is empty, or
/// [`CapacitorError::InvalidParameter`] if any entry is not strictly
/// positive and finite.
///
/// # Examples
///
/// ```
/// use valenx_capacitor::network::series;
///
/// // Two equal 2 uF capacitors in series give 1 uF.
/// let c = series(&[2.0e-6, 2.0e-6]).unwrap();
/// assert!((c - 1.0e-6).abs() < 1e-18);
/// ```
pub fn series(caps_f: &[f64]) -> Result<f64> {
    if caps_f.is_empty() {
        return Err(CapacitorError::EmptyNetwork(
            "series combination needs at least one capacitor",
        ));
    }
    let mut recip_sum = 0.0;
    for (i, &c) in caps_f.iter().enumerate() {
        if !(c.is_finite() && c > 0.0) {
            return Err(CapacitorError::InvalidParameter {
                name: index_name(i),
                value: c,
                reason: "each branch capacitance must be positive",
            });
        }
        recip_sum += 1.0 / c;
    }
    Ok(1.0 / recip_sum)
}

/// Map a branch index to one of a small set of `'static` names for error
/// reporting, so [`CapacitorError::InvalidParameter`] can keep a
/// `&'static str` field without allocating per call.
fn index_name(i: usize) -> &'static str {
    const NAMES: [&str; 8] = [
        "caps_f[0]",
        "caps_f[1]",
        "caps_f[2]",
        "caps_f[3]",
        "caps_f[4]",
        "caps_f[5]",
        "caps_f[6]",
        "caps_f[7]",
    ];
    NAMES.get(i).copied().unwrap_or("caps_f[..]")
}
