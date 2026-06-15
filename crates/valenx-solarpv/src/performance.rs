//! Operating-point metrics: maximum power point, fill factor, efficiency.
//!
//! These functions take a constructed [`SingleDiode`] and reduce its
//! I-V curve to the figures of merit quoted on a PV datasheet:
//!
//! - the **maximum power point** (MPP): the `(Vmp, Imp, Pmax)` that
//!   maximises `P(V) = V * I(V)` over `0 <= V <= Voc`;
//! - the **fill factor** `FF = Pmax / (Voc * Isc)`, a dimensionless
//!   curve-squareness measure that for good cells sits in roughly
//!   `0.7 - 0.85`;
//! - the **module efficiency** `eta = Pmax / (irradiance * area)`, the
//!   fraction of incident optical power converted to electrical power.

use crate::diode::SingleDiode;
use crate::error::{Result, SolarPvError};

/// The maximum-power operating point of a cell.
///
/// All fields are at the point that maximises `V * I` on the
/// power-producing branch `0 <= V <= Voc`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct MaxPowerPoint {
    /// Voltage at maximum power `Vmp`, in volts.
    pub v_mp: f64,
    /// Current at maximum power `Imp`, in amperes.
    pub i_mp: f64,
    /// Maximum power `Pmax = Vmp * Imp`, in watts.
    pub p_max: f64,
}

/// Locate the maximum power point of `cell` by a coarse scan over
/// `[0, Voc]` followed by a golden-section refinement of the bracket
/// around the best scan sample.
///
/// `samples` controls the coarse scan resolution; values around
/// `200`-`1000` give a tight bracket cheaply. The refinement then drives
/// the voltage bracket below a 1e-9 V width, so the reported `Pmax` is
/// accurate to well past datasheet precision.
///
/// # Errors
///
/// Returns [`SolarPvError::Invalid`] if `samples < 2`. Propagates any
/// error from [`SingleDiode::voc`] or [`SingleDiode::current_at`]
/// (e.g. a non-convergent root-find, or a zero-saturation-current cell
/// whose `Voc` is undefined).
pub fn max_power_point(cell: &SingleDiode, samples: usize) -> Result<MaxPowerPoint> {
    if samples < 2 {
        return Err(SolarPvError::invalid(
            "samples",
            format!("need at least 2 scan samples, got {samples}"),
        ));
    }
    let voc = cell.voc()?;
    if voc <= 0.0 {
        // Dark / degenerate cell: the only operating point is the origin.
        return Ok(MaxPowerPoint {
            v_mp: 0.0,
            i_mp: cell.current_at(0.0)?,
            p_max: 0.0,
        });
    }

    // Coarse scan: find the sample voltage with the largest power, and
    // keep its immediate neighbours to bracket the true peak.
    let mut best_k = 0usize;
    let mut best_p = f64::NEG_INFINITY;
    for k in 0..=samples {
        let v = voc * (k as f64) / (samples as f64);
        let p = cell.power_at(v)?;
        if p > best_p {
            best_p = p;
            best_k = k;
        }
    }
    let dv = voc / (samples as f64);
    let mut lo = voc.min((best_k as f64).max(1.0) * dv - dv).max(0.0);
    let mut hi = ((best_k as f64) * dv + dv).min(voc);
    if hi <= lo {
        // Peak sits at an endpoint; widen minimally to a valid bracket.
        lo = (lo - dv).max(0.0);
        hi = (hi + dv).min(voc);
    }

    // Golden-section search for the maximiser of power on [lo, hi].
    const INV_PHI: f64 = 0.618_033_988_749_894_8; // 1/phi
    let mut a = lo;
    let mut b = hi;
    let mut c = b - (b - a) * INV_PHI;
    let mut d = a + (b - a) * INV_PHI;
    let mut fc = cell.power_at(c)?;
    let mut fd = cell.power_at(d)?;
    const MAX_ITERS: u32 = 200;
    const TOL: f64 = 1e-9;
    for _ in 0..MAX_ITERS {
        if (b - a).abs() <= TOL {
            break;
        }
        if fc > fd {
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) * INV_PHI;
            fc = cell.power_at(c)?;
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) * INV_PHI;
            fd = cell.power_at(d)?;
        }
    }
    let v_mp = 0.5 * (a + b);
    let i_mp = cell.current_at(v_mp)?;
    let p_max = v_mp * i_mp;
    Ok(MaxPowerPoint { v_mp, i_mp, p_max })
}

/// Fill factor `FF = Pmax / (Voc * Isc)` (dimensionless).
///
/// Measures how "square" the I-V curve is: the ratio of the largest
/// rectangle that fits under the curve (`Pmax`) to the bounding
/// rectangle (`Voc * Isc`). For a physical cell the maximum power is
/// always strictly inside the bounding box, so `FF` lies in `(0, 1)`;
/// good crystalline-silicon cells sit around `0.7 - 0.85`.
///
/// # Errors
///
/// Returns [`SolarPvError::Undefined`] if `Voc * Isc == 0` (a dark or
/// degenerate cell has no defined fill factor). Propagates errors from
/// the underlying [`max_power_point`], [`SingleDiode::voc`] and
/// [`SingleDiode::isc`].
pub fn fill_factor(cell: &SingleDiode, samples: usize) -> Result<f64> {
    let voc = cell.voc()?;
    let isc = cell.isc()?;
    let denom = voc * isc;
    if denom == 0.0 {
        return Err(SolarPvError::undefined(
            "fill_factor",
            "Voc * Isc is zero (dark or degenerate cell)",
        ));
    }
    let mpp = max_power_point(cell, samples)?;
    Ok(mpp.p_max / denom)
}

/// Module power-conversion efficiency `eta = Pmax / (irradiance * area)`
/// (dimensionless fraction, multiply by 100 for a percentage).
///
/// `irradiance_w_per_m2` is the plane-of-array irradiance in W/m^2
/// (1000 at "one sun" / STC) and `area_m2` the active cell/module area
/// in square metres. The denominator `irradiance * area` is the incident
/// optical power in watts; the result is the fraction of that power
/// delivered at the maximum power point.
///
/// # Errors
///
/// Returns [`SolarPvError::Invalid`] if `irradiance_w_per_m2 <= 0` or
/// `area_m2 <= 0` (non-finite included), and [`SolarPvError::Undefined`]
/// only via the propagated MPP path. Propagates any error from
/// [`max_power_point`].
pub fn efficiency(
    cell: &SingleDiode,
    irradiance_w_per_m2: f64,
    area_m2: f64,
    samples: usize,
) -> Result<f64> {
    if !irradiance_w_per_m2.is_finite() || irradiance_w_per_m2 <= 0.0 {
        return Err(SolarPvError::invalid(
            "irradiance_w_per_m2",
            format!("irradiance must be finite and > 0, got {irradiance_w_per_m2}"),
        ));
    }
    if !area_m2.is_finite() || area_m2 <= 0.0 {
        return Err(SolarPvError::invalid(
            "area_m2",
            format!("area must be finite and > 0, got {area_m2}"),
        ));
    }
    let incident_power_w = irradiance_w_per_m2 * area_m2;
    let mpp = max_power_point(cell, samples)?;
    Ok(mpp.p_max / incident_power_w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{STC_IRRADIANCE_W_PER_M2, STC_TEMPERATURE_K};

    fn si_cell() -> SingleDiode {
        SingleDiode::new(3.8, 1.0e-9, 1.2, STC_TEMPERATURE_K, 0.01, 200.0).unwrap()
    }

    fn ideal_cell() -> SingleDiode {
        SingleDiode::ideal(3.8, 1.0e-9, 1.2, STC_TEMPERATURE_K).unwrap()
    }

    /// VALIDATE: Pmax <= Voc * Isc (the MPP rectangle fits in the
    /// bounding box), and the MPP voltage/current sit inside the curve.
    #[test]
    fn pmax_does_not_exceed_voc_times_isc() {
        for c in [ideal_cell(), si_cell()] {
            let voc = c.voc().unwrap();
            let isc = c.isc().unwrap();
            let mpp = max_power_point(&c, 400).unwrap();
            assert!(
                mpp.p_max <= voc * isc + 1e-9,
                "Pmax={} > Voc*Isc={}",
                mpp.p_max,
                voc * isc
            );
            assert!(mpp.v_mp > 0.0 && mpp.v_mp < voc, "Vmp = {}", mpp.v_mp);
            assert!(
                mpp.i_mp > 0.0 && mpp.i_mp < isc + 1e-9,
                "Imp = {}",
                mpp.i_mp
            );
        }
    }

    /// The located MPP is genuinely a maximum: power at Vmp is at least as
    /// large as power at a dense set of other voltages on [0, Voc].
    #[test]
    fn mpp_is_a_true_maximum() {
        let c = si_cell();
        let voc = c.voc().unwrap();
        let mpp = max_power_point(&c, 500).unwrap();
        for k in 0..=200 {
            let v = voc * (k as f64) / 200.0;
            let p = c.power_at(v).unwrap();
            assert!(
                p <= mpp.p_max + 1e-6,
                "found higher power {p} at v={v} than Pmax={}",
                mpp.p_max
            );
        }
    }

    /// Pmax reported equals Vmp * Imp recomputed from the diode model.
    #[test]
    fn pmax_is_consistent_with_vmp_imp() {
        let c = si_cell();
        let mpp = max_power_point(&c, 400).unwrap();
        let i_check = c.current_at(mpp.v_mp).unwrap();
        assert!((mpp.i_mp - i_check).abs() < 1e-9, "Imp inconsistent");
        assert!(
            (mpp.p_max - mpp.v_mp * mpp.i_mp).abs() < 1e-12,
            "Pmax != Vmp*Imp"
        );
    }

    /// VALIDATE: FF in (0, 1), and for these realistic parameters in the
    /// physically expected 0.7-0.85 band.
    #[test]
    fn fill_factor_in_typical_band() {
        let c = si_cell();
        let ff = fill_factor(&c, 600).unwrap();
        assert!(ff > 0.0 && ff < 1.0, "FF = {ff}");
        assert!((0.70..=0.86).contains(&ff), "FF = {ff} outside 0.70-0.86");
    }

    /// FF equals Pmax / (Voc * Isc) computed independently from the
    /// component metrics.
    #[test]
    fn fill_factor_matches_definition() {
        let c = si_cell();
        let voc = c.voc().unwrap();
        let isc = c.isc().unwrap();
        let mpp = max_power_point(&c, 600).unwrap();
        let ff_direct = mpp.p_max / (voc * isc);
        let ff = fill_factor(&c, 600).unwrap();
        assert!((ff - ff_direct).abs() < 1e-9, "ff={ff}, direct={ff_direct}");
    }

    /// VALIDATE the efficiency formula: eta = Pmax / (G * A) recomputed
    /// from the MPP, and lands in a sane single-junction range.
    #[test]
    fn efficiency_matches_formula() {
        let c = si_cell();
        // Pick an area so Pmax/(G*A) is a believable cell efficiency.
        let area = 0.0243; // m^2, ~156mm pseudo-square wafer
        let mpp = max_power_point(&c, 600).unwrap();
        let expected = mpp.p_max / (STC_IRRADIANCE_W_PER_M2 * area);
        let eta = efficiency(&c, STC_IRRADIANCE_W_PER_M2, area, 600).unwrap();
        assert!((eta - expected).abs() < 1e-12, "eta = {eta}");
        assert!(eta > 0.0 && eta < 1.0, "eta = {eta}");
    }

    /// Efficiency scales inversely with area for fixed irradiance, and
    /// inversely with irradiance for fixed area (Pmax fixed by the cell).
    #[test]
    fn efficiency_scales_with_denominator() {
        let c = si_cell();
        let e1 = efficiency(&c, STC_IRRADIANCE_W_PER_M2, 0.02, 600).unwrap();
        let e2 = efficiency(&c, STC_IRRADIANCE_W_PER_M2, 0.04, 600).unwrap();
        // Double the area -> half the efficiency (same Pmax, since the
        // cell's Iph is fixed here and independent of the quoted area).
        assert!((e2 - 0.5 * e1).abs() < 1e-9, "e1={e1}, e2={e2}");
    }

    #[test]
    fn efficiency_rejects_non_positive_inputs() {
        let c = si_cell();
        assert!(efficiency(&c, 0.0, 1.0, 100).is_err());
        assert!(efficiency(&c, -100.0, 1.0, 100).is_err());
        assert!(efficiency(&c, 1000.0, 0.0, 100).is_err());
        assert!(efficiency(&c, 1000.0, -1.0, 100).is_err());
        assert!(efficiency(&c, f64::NAN, 1.0, 100).is_err());
    }

    #[test]
    fn max_power_point_rejects_too_few_samples() {
        let c = si_cell();
        assert!(max_power_point(&c, 1).is_err());
        assert!(max_power_point(&c, 0).is_err());
    }

    #[test]
    fn fill_factor_undefined_for_dark_cell() {
        let dark = SingleDiode::ideal(0.0, 1e-9, 1.0, 300.0).unwrap();
        // Voc = 0 -> Voc*Isc = 0 -> undefined.
        assert!(fill_factor(&dark, 100).is_err());
    }
}
