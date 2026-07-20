//! The area-Mach relation for a converging-diverging (de Laval) nozzle and
//! its sub/supersonic inversion.
//!
//! For steady isentropic flow of a perfect gas, continuity plus the
//! stagnation relations tie the local cross-sectional area `A` to the
//! sonic-throat area `A*` purely through the local Mach number `M`:
//!
//! ```text
//! A/A* = (1/M) * [ (2/(g+1)) * (1 + (g-1)/2 * M^2) ] ^ ( (g+1) / (2*(g-1)) )
//! ```
//!
//! where `g` is `gamma`. The function has a single global minimum of `1.0`
//! at `M = 1` (the throat is sonic when the nozzle is choked) and rises
//! toward `+inf` as `M -> 0` and as `M -> inf`. Consequently every area
//! ratio `A/A* > 1` corresponds to **two** Mach numbers — one subsonic
//! (`M < 1`, in the converging section) and one supersonic (`M > 1`, in the
//! diverging section). [`mach_from_area_ratio`] selects the branch.

use crate::error::{check_gamma, check_mach, NozzleError, Result};

/// The area-Mach ratio `A/A*` at a given Mach number.
///
/// See the [module docs](crate::area) for the governing equation. At
/// `M = 1` this returns exactly `1.0` (the choked throat); it diverges as
/// `M -> 0`.
///
/// # Errors
///
/// Returns a [`NozzleError`] for an out-of-range `gamma`, a negative /
/// non-finite `mach`, or `mach == 0` (where `A/A*` is infinite and so is
/// rejected as out of domain).
///
/// # Example
///
/// ```
/// // The defining ground truth: M = 1 gives A/A* = 1.
/// let r = valenx_nozzle::area_ratio(1.0, 1.4).unwrap();
/// assert!((r - 1.0).abs() < 1e-12);
/// ```
pub fn area_ratio(mach: f64, gamma: f64) -> Result<f64> {
    check_gamma(gamma)?;
    check_mach(mach)?;
    if mach == 0.0 {
        // A/A* -> +inf at M = 0; treat as out of domain rather than return inf.
        return Err(NozzleError::NegativeMach { value: mach });
    }
    let g = gamma;
    let term = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * mach * mach);
    let exponent = (g + 1.0) / (2.0 * (g - 1.0));
    Ok(term.powf(exponent) / mach)
}

/// Which root of the area-Mach relation to return.
///
/// For any `A/A* > 1` there is a subsonic and a supersonic Mach number with
/// the same area ratio; this picks the physical branch for the section of
/// the nozzle you are in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Branch {
    /// The converging-section root, `0 < M < 1`.
    Subsonic,
    /// The diverging-section root, `M > 1`.
    Supersonic,
}

/// Recover the Mach number from an area ratio `A/A*` on a chosen branch.
///
/// Solves the area-Mach relation for `M` by Newton's method. At
/// `A/A* == 1` both branches return `M = 1` (the choked throat). The
/// iteration is seeded from a closed-form approximation so it converges in a
/// handful of steps over the usual nozzle range.
///
/// # Errors
///
/// Returns [`NozzleError::AreaRatioBelowChoke`] if `area_ratio < 1`,
/// [`NozzleError::NotFinite`] if it is `NaN`/`±∞`, a validation error for a
/// bad `gamma`, or [`NozzleError::NotConverged`] if Newton fails to settle.
///
/// # Example
///
/// ```
/// use valenx_nozzle::{mach_from_area_ratio, Branch};
/// // Air, A/A* = 2: the supersonic root is ~2.1972, the subsonic ~0.3059.
/// let sup = mach_from_area_ratio(2.0, 1.4, Branch::Supersonic).unwrap();
/// let sub = mach_from_area_ratio(2.0, 1.4, Branch::Subsonic).unwrap();
/// assert!((sup - 2.197_198_1).abs() < 1e-6);
/// assert!((sub - 0.305_903_8).abs() < 1e-6);
/// ```
pub fn mach_from_area_ratio(area_ratio: f64, gamma: f64, branch: Branch) -> Result<f64> {
    check_gamma(gamma)?;
    if !area_ratio.is_finite() {
        return Err(NozzleError::NotFinite {
            name: "area_ratio",
            value: area_ratio,
        });
    }
    if area_ratio < 1.0 {
        return Err(NozzleError::AreaRatioBelowChoke { value: area_ratio });
    }
    // The throat: both branches meet at M = 1. Handle exactly so the
    // 1/(M-derivative) Newton step (which is flat here) is never taken.
    if (area_ratio - 1.0).abs() < 1e-12 {
        return Ok(1.0);
    }

    let g = gamma;
    let exponent = (g + 1.0) / (2.0 * (g - 1.0));

    // f(M) = (2/(g+1) * (1 + (g-1)/2 M^2))^exp / M  -  area_ratio.
    let f = |m: f64| -> f64 {
        let term = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * m * m);
        term.powf(exponent) / m - area_ratio
    };

    // Initial guess. These are the standard subsonic/supersonic seeds; the
    // exact constant doesn't matter much because Newton converges fast, but
    // good seeds keep us inside the monotone region of each branch.
    let mut m = match branch {
        Branch::Subsonic => {
            // A smooth subsonic approximation valid for moderate area
            // ratios; bounded into (0, 1).
            let guess = 1.0 / (area_ratio + 0.5);
            guess.clamp(1e-4, 0.999_9)
        }
        Branch::Supersonic => {
            // Grows with area ratio; kept above 1.
            (1.0 + (area_ratio - 1.0)).max(1.000_1)
        }
    };

    // Newton with a numerical derivative; the analytic derivative is messy
    // and a central difference is plenty here.
    let max_iters: u32 = 100;
    for i in 0..max_iters {
        let fm = f(m);
        if fm.abs() < 1e-13 {
            return Ok(m);
        }
        let h = 1e-7 * m.max(1.0);
        let dfm = (f(m + h) - f(m - h)) / (2.0 * h);
        if dfm == 0.0 || !dfm.is_finite() {
            return Err(NozzleError::NotConverged {
                residual: fm.abs(),
                iters: i,
            });
        }
        let mut next = m - fm / dfm;
        // Keep each branch on its own side of the sonic point so a long
        // Newton step can't jump across M = 1 into the other root.
        match branch {
            Branch::Subsonic => next = next.clamp(1e-6, 0.999_999),
            Branch::Supersonic => {
                if next <= 1.0 {
                    next = 1.0 + (m - 1.0) * 0.5;
                }
            }
        }
        if (next - m).abs() < 1e-12 {
            m = next;
            break;
        }
        m = next;
    }

    let residual = f(m).abs();
    if residual < 1e-9 {
        Ok(m)
    } else {
        Err(NozzleError::NotConverged {
            residual,
            iters: max_iters,
        })
    }
}

/// The dimensionless **mass-flow (flow-rate) parameter** at a given Mach.
///
/// This is the grouping `m_dot * sqrt(T0) / (A * p0)` divided by
/// `sqrt(gamma / R)` — i.e. the part that depends only on `M` and `gamma`:
///
/// ```text
/// MFP(M) = M * sqrt(g) * (1 + (g-1)/2 M^2) ^ ( -(g+1) / (2*(g-1)) )
/// ```
///
/// It rises from `0` at `M = 0` to its maximum at the sonic point `M = 1`
/// (the choked value), then falls again — the analytic expression of the
/// fact that a choked throat passes the maximum possible mass flow for a
/// given stagnation state and area.
///
/// # Errors
///
/// Returns a [`NozzleError`] for an out-of-range `gamma` or a negative /
/// non-finite `mach`.
///
/// # Example
///
/// ```
/// // The parameter peaks at M = 1 and is symmetric-ish in value either side.
/// let at_throat = valenx_nozzle::mass_flow_parameter(1.0, 1.4).unwrap();
/// let subsonic = valenx_nozzle::mass_flow_parameter(0.8, 1.4).unwrap();
/// let supersonic = valenx_nozzle::mass_flow_parameter(1.2, 1.4).unwrap();
/// assert!(at_throat > subsonic && at_throat > supersonic);
/// ```
pub fn mass_flow_parameter(mach: f64, gamma: f64) -> Result<f64> {
    check_gamma(gamma)?;
    check_mach(mach)?;
    let g = gamma;
    let base = 1.0 + 0.5 * (g - 1.0) * mach * mach;
    let exponent = -(g + 1.0) / (2.0 * (g - 1.0));
    Ok(mach * g.sqrt() * base.powf(exponent))
}
