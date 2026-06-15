//! Heat-pump *balance point*: the outdoor temperature at which the
//! building's heating load exactly equals the heat pump's heating
//! capacity.
//!
//! ## The textbook picture
//!
//! For an air-source heat pump in winter, two straight lines on a
//! (outdoor-temperature, kilowatts) plot tell the story.
//!
//! The building **load** rises as it gets colder. A building loses heat
//! in proportion to the indoor-minus-outdoor temperature difference, so
//! the load is a line that is zero at the balanced indoor set-point
//! `T_balance` and increases with slope `ua` (the building's overall loss
//! coefficient, kW per kelvin) as the outdoor temperature `T` drops below
//! it:
//!
//! ```text
//! load(T) = ua * (T_balance - T)        for T < T_balance
//! ```
//!
//! The heat-pump **capacity** *falls* as it gets colder. An air-source
//! unit pulls heat from outdoor air, and there is less of it to pull when
//! the air is cold, so the rated capacity drops roughly linearly:
//!
//! ```text
//! capacity(T) = cap_ref + slope * (T - T_ref)
//! ```
//!
//! with a *positive* `slope` (capacity goes down as `T` goes down).
//!
//! The **balance point** `T_bp` is where the two lines cross,
//! `load(T_bp) = capacity(T_bp)`. Above it the heat pump has spare
//! capacity; below it the load wins and a backup (resistive / gas) heat
//! source must make up the deficit. Because both curves are affine, the
//! crossing has a closed form — but this module *also* finds it with a
//! bracketed bisection root-find on the residual `capacity - load`, so
//! the same machinery extends to the non-linear capacity tables real
//! tools use. The closed form is kept as a cross-check.

use serde::{Deserialize, Serialize};

use crate::error::{HeatPumpError, Result};

/// A building heating-load line: zero at `t_balance_c`, rising with
/// slope `ua_kw_per_k` as the outdoor temperature falls below it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct LoadLine {
    /// Indoor balance temperature (degrees Celsius) at and above which
    /// the heating load is zero — free gains cover the loss.
    pub t_balance_c: f64,
    /// Overall building heat-loss coefficient `UA`, in kilowatts per
    /// kelvin. Must be strictly positive.
    pub ua_kw_per_k: f64,
}

impl LoadLine {
    /// Construct a validated load line.
    ///
    /// # Errors
    ///
    /// [`HeatPumpError::Invalid`] if `t_balance_c` is not finite, or if
    /// `ua_kw_per_k` is not finite or not strictly positive.
    pub fn new(t_balance_c: f64, ua_kw_per_k: f64) -> Result<Self> {
        if !t_balance_c.is_finite() {
            return Err(HeatPumpError::invalid(
                "t_balance_c",
                "must be a finite temperature",
            ));
        }
        if !ua_kw_per_k.is_finite() || ua_kw_per_k <= 0.0 {
            return Err(HeatPumpError::invalid(
                "ua_kw_per_k",
                format!("must be a finite positive loss coefficient, got {ua_kw_per_k}"),
            ));
        }
        Ok(Self {
            t_balance_c,
            ua_kw_per_k,
        })
    }

    /// The heating load, in kilowatts, at outdoor temperature
    /// `t_outdoor_c` (degrees Celsius).
    ///
    /// Returns `0` at or above [`t_balance_c`](Self::t_balance_c) and a
    /// positive, increasing value below it.
    pub fn load_kw(&self, t_outdoor_c: f64) -> f64 {
        (self.t_balance_c - t_outdoor_c).max(0.0) * self.ua_kw_per_k
    }
}

/// An air-source heat pump's heating-capacity line, anchored at a
/// reference rating point and *falling* with a positive `slope` as the
/// outdoor temperature drops.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CapacityLine {
    /// Reference outdoor temperature (degrees Celsius) at which the unit
    /// is rated.
    pub t_ref_c: f64,
    /// Rated heating capacity (kilowatts) at [`t_ref_c`](Self::t_ref_c).
    /// Must be strictly positive.
    pub cap_ref_kw: f64,
    /// Capacity slope, in kilowatts per kelvin. Stored as the magnitude
    /// of the *fall* per kelvin of cooling: capacity = `cap_ref + slope *
    /// (T - T_ref)`, so a positive slope means capacity drops as `T`
    /// drops. Must be strictly positive (a real air-source unit always
    /// loses capacity in the cold).
    pub slope_kw_per_k: f64,
}

impl CapacityLine {
    /// Construct a validated capacity line.
    ///
    /// # Errors
    ///
    /// [`HeatPumpError::Invalid`] if any field is not finite, if
    /// `cap_ref_kw` is not strictly positive, or if `slope_kw_per_k` is
    /// not strictly positive.
    pub fn new(t_ref_c: f64, cap_ref_kw: f64, slope_kw_per_k: f64) -> Result<Self> {
        if !t_ref_c.is_finite() {
            return Err(HeatPumpError::invalid(
                "t_ref_c",
                "must be a finite temperature",
            ));
        }
        if !cap_ref_kw.is_finite() || cap_ref_kw <= 0.0 {
            return Err(HeatPumpError::invalid(
                "cap_ref_kw",
                format!("must be a finite positive capacity, got {cap_ref_kw}"),
            ));
        }
        if !slope_kw_per_k.is_finite() || slope_kw_per_k <= 0.0 {
            return Err(HeatPumpError::invalid(
                "slope_kw_per_k",
                format!("must be a finite positive slope, got {slope_kw_per_k}"),
            ));
        }
        Ok(Self {
            t_ref_c,
            cap_ref_kw,
            slope_kw_per_k,
        })
    }

    /// The available heating capacity, in kilowatts, at outdoor
    /// temperature `t_outdoor_c` (degrees Celsius).
    ///
    /// Clamped to be non-negative — below the temperature at which the
    /// extrapolated line would go negative the unit simply delivers
    /// nothing.
    pub fn capacity_kw(&self, t_outdoor_c: f64) -> f64 {
        (self.cap_ref_kw + self.slope_kw_per_k * (t_outdoor_c - self.t_ref_c)).max(0.0)
    }
}

/// The result of a balance-point solve: the crossover temperature and
/// the (equal) load and capacity there.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BalancePoint {
    /// Outdoor temperature (degrees Celsius) at which load equals
    /// capacity.
    pub t_balance_c: f64,
    /// The heating load (kilowatts) at the balance point — equal, to
    /// solver tolerance, to [`capacity_kw`](Self::capacity_kw).
    pub load_kw: f64,
    /// The heating capacity (kilowatts) at the balance point.
    pub capacity_kw: f64,
}

impl BalancePoint {
    /// The signed residual `capacity - load` at the solved point — should
    /// be within the solver tolerance of zero.
    pub fn residual_kw(&self) -> f64 {
        self.capacity_kw - self.load_kw
    }
}

/// Find the balance point by bracketed bisection on the residual
/// `capacity(T) - load(T)` over the outdoor-temperature interval
/// `[t_lo_c, t_hi_c]`.
///
/// The residual is *decreasing* in `T` only in the usual regime, but the
/// solver makes no monotonicity assumption: it requires only a sign
/// change between the endpoints and bisects to `tol_kw` on the residual
/// (or until the bracket is narrower than `1e-12` °C). This is the same
/// root-find that drives the [`balance_point_linear`] closed form below,
/// kept general so it survives a swap to a non-linear capacity table.
///
/// # Errors
///
/// Returns [`HeatPumpError::Invalid`] if the interval is empty / not
/// finite or `tol_kw` is not finite and positive. Returns
/// [`HeatPumpError::NoConvergence`] if the residual does not change sign
/// across `[t_lo_c, t_hi_c]` (the lines never cross inside the bracket —
/// capacity dominates or load dominates throughout).
///
/// # Examples
///
/// ```
/// use valenx_heatpump::balance::{solve_balance_point, CapacityLine, LoadLine};
/// let load = LoadLine::new(18.0, 0.5).unwrap();
/// let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();
/// let bp = solve_balance_point(&load, &cap, -25.0, 18.0, 1e-9).unwrap();
/// assert!(bp.residual_kw().abs() < 1e-6);
/// ```
pub fn solve_balance_point(
    load: &LoadLine,
    capacity: &CapacityLine,
    t_lo_c: f64,
    t_hi_c: f64,
    tol_kw: f64,
) -> Result<BalancePoint> {
    if !t_lo_c.is_finite() || !t_hi_c.is_finite() || t_lo_c >= t_hi_c {
        return Err(HeatPumpError::invalid(
            "interval",
            format!("require finite t_lo_c < t_hi_c, got [{t_lo_c}, {t_hi_c}]"),
        ));
    }
    if !tol_kw.is_finite() || tol_kw <= 0.0 {
        return Err(HeatPumpError::invalid(
            "tol_kw",
            format!("must be a finite positive tolerance, got {tol_kw}"),
        ));
    }

    let residual = |t: f64| capacity.capacity_kw(t) - load.load_kw(t);
    let mut lo = t_lo_c;
    let mut hi = t_hi_c;
    let mut f_lo = residual(lo);
    let f_hi = residual(hi);

    // Exact hit at an endpoint.
    if f_lo.abs() <= tol_kw {
        return Ok(point_at(load, capacity, lo));
    }
    if f_hi.abs() <= tol_kw {
        return Ok(point_at(load, capacity, hi));
    }
    if f_lo.signum() == f_hi.signum() {
        return Err(HeatPumpError::no_convergence(format!(
            "residual capacity-load does not change sign across [{t_lo_c}, {t_hi_c}] \
             (endpoints {f_lo} and {f_hi} kW)"
        )));
    }

    // Bisection.
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let f_mid = residual(mid);
        if f_mid.abs() <= tol_kw || (hi - lo) < 1e-12 {
            return Ok(point_at(load, capacity, mid));
        }
        if f_mid.signum() == f_lo.signum() {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
        }
    }
    Ok(point_at(load, capacity, 0.5 * (lo + hi)))
}

/// Build a [`BalancePoint`] by sampling the load and capacity at a solved
/// temperature.
fn point_at(load: &LoadLine, capacity: &CapacityLine, t_c: f64) -> BalancePoint {
    BalancePoint {
        t_balance_c: t_c,
        load_kw: load.load_kw(t_c),
        capacity_kw: capacity.capacity_kw(t_c),
    }
}

/// The closed-form balance point for the affine load and capacity lines,
/// ignoring the non-negativity clamps.
///
/// Setting `ua * (T_balance - T) = cap_ref + slope * (T - T_ref)` and
/// solving for `T` gives:
///
/// ```text
///        ua * T_balance - cap_ref + slope * T_ref
/// T_bp = ----------------------------------------
///                     ua + slope
/// ```
///
/// This is the analytic cross-check for [`solve_balance_point`]; the
/// denominator `ua + slope` is strictly positive because both terms are
/// validated positive, so the division is always well defined.
///
/// # Errors
///
/// This never fails for validated [`LoadLine`] / [`CapacityLine`] inputs
/// (the denominator cannot be zero), but returns a [`Result`] for a
/// uniform signature.
///
/// # Examples
///
/// ```
/// use valenx_heatpump::balance::{balance_point_linear, CapacityLine, LoadLine};
/// let load = LoadLine::new(18.0, 0.5).unwrap();
/// let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();
/// let t = balance_point_linear(&load, &cap).unwrap();
/// assert!(t.is_finite());
/// ```
pub fn balance_point_linear(load: &LoadLine, capacity: &CapacityLine) -> Result<f64> {
    let denom = load.ua_kw_per_k + capacity.slope_kw_per_k;
    // Both terms are validated strictly positive, so `denom > 0`; the
    // guard is defensive only.
    if denom <= 0.0 {
        return Err(HeatPumpError::no_convergence(
            "load slope + capacity slope is non-positive",
        ));
    }
    let numer = load.ua_kw_per_k * load.t_balance_c - capacity.cap_ref_kw
        + capacity.slope_kw_per_k * capacity.t_ref_c;
    Ok(numer / denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for floating-point comparisons.
    const EPS: f64 = 1e-9;

    /// Ground truth, fully hand-computed.
    ///
    /// Load: `T_balance = 18 °C`, `UA = 0.5 kW/K`.
    /// Capacity: rated `10 kW` at `8.3 °C`, falling `0.25 kW/K`.
    ///
    /// Closed form:
    /// `T_bp = (0.5*18 - 10 + 0.25*8.3) / (0.5 + 0.25)`
    ///       = `(9 - 10 + 2.075) / 0.75`
    ///       = `1.075 / 0.75`
    ///       = `1.4333333... °C`.
    ///
    /// At that T: `load = 0.5*(18 - 1.43333) = 8.28333 kW`,
    /// and `capacity = 10 + 0.25*(1.43333 - 8.3) = 8.28333 kW` — equal.
    #[test]
    fn balance_point_matches_hand_computed_ground_truth() {
        let load = LoadLine::new(18.0, 0.5).unwrap();
        let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();

        let expected_t = 1.075 / 0.75;
        let t_closed = balance_point_linear(&load, &cap).unwrap();
        assert!(
            (t_closed - expected_t).abs() < EPS,
            "closed form: {t_closed}"
        );

        let bp = solve_balance_point(&load, &cap, -25.0, 18.0, 1e-12).unwrap();
        assert!((bp.t_balance_c - expected_t).abs() < 1e-6);
        let expected_kw = 0.5 * (18.0 - expected_t);
        assert!((bp.load_kw - expected_kw).abs() < 1e-6);
        assert!((bp.capacity_kw - expected_kw).abs() < 1e-6);
    }

    /// VALIDATE: at the solved balance point load equals capacity.
    #[test]
    fn load_equals_capacity_at_balance_point() {
        let load = LoadLine::new(20.0, 0.42).unwrap();
        let cap = CapacityLine::new(7.0, 12.0, 0.3).unwrap();
        let bp = solve_balance_point(&load, &cap, -30.0, 20.0, 1e-10).unwrap();
        assert!(
            (bp.load_kw - bp.capacity_kw).abs() < 1e-6,
            "load {} != capacity {}",
            bp.load_kw,
            bp.capacity_kw
        );
        assert!(bp.residual_kw().abs() < 1e-6);
    }

    /// The bisection solver and the closed form must agree across a sweep
    /// of building / unit combinations.
    #[test]
    fn bisection_agrees_with_closed_form() {
        let cases = [
            (18.0, 0.5, 8.3, 10.0, 0.25),
            (21.0, 0.30, 7.0, 14.0, 0.40),
            (16.0, 0.80, 2.0, 9.0, 0.10),
            (22.0, 0.25, -5.0, 18.0, 0.5),
        ];
        for (tb, ua, tref, capref, slope) in cases {
            let load = LoadLine::new(tb, ua).unwrap();
            let cap = CapacityLine::new(tref, capref, slope).unwrap();
            let t_closed = balance_point_linear(&load, &cap).unwrap();
            let bp = solve_balance_point(&load, &cap, -40.0, tb, 1e-12).unwrap();
            assert!(
                (bp.t_balance_c - t_closed).abs() < 1e-5,
                "bisection {} vs closed {t_closed} for case ua={ua}",
                bp.t_balance_c
            );
        }
    }

    /// Above the balance point the heat pump has spare capacity
    /// (capacity > load); below it the load wins (capacity < load).
    #[test]
    fn capacity_dominates_above_and_load_dominates_below() {
        let load = LoadLine::new(18.0, 0.5).unwrap();
        let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();
        let t_bp = balance_point_linear(&load, &cap).unwrap();

        let above = t_bp + 5.0;
        assert!(cap.capacity_kw(above) > load.load_kw(above));

        let below = t_bp - 5.0;
        assert!(cap.capacity_kw(below) < load.load_kw(below));
    }

    /// Load rises as the outdoor temperature falls; capacity falls.
    #[test]
    fn load_rises_and_capacity_falls_as_it_gets_colder() {
        let load = LoadLine::new(18.0, 0.5).unwrap();
        let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();
        assert!(load.load_kw(-10.0) > load.load_kw(5.0));
        assert!(cap.capacity_kw(-10.0) < cap.capacity_kw(5.0));
        // No load demanded at/above the balance temperature.
        assert!((load.load_kw(18.0)).abs() < EPS);
        assert!((load.load_kw(25.0)).abs() < EPS);
    }

    /// If the lines never cross inside the bracket the solver reports a
    /// non-convergence rather than returning a bogus root.
    #[test]
    fn non_crossing_bracket_reports_no_convergence() {
        // Huge unit, tiny load: capacity dominates the whole interval.
        let load = LoadLine::new(18.0, 0.01).unwrap();
        let cap = CapacityLine::new(8.3, 100.0, 0.1).unwrap();
        let err = solve_balance_point(&load, &cap, 0.0, 18.0, 1e-9).unwrap_err();
        assert!(matches!(err, HeatPumpError::NoConvergence { .. }));
        assert_eq!(err.code(), "heatpump.no_convergence");
    }

    #[test]
    fn bad_interval_and_tolerance_are_rejected() {
        let load = LoadLine::new(18.0, 0.5).unwrap();
        let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();
        // Inverted interval.
        assert!(solve_balance_point(&load, &cap, 18.0, -25.0, 1e-9).is_err());
        // Non-positive tolerance.
        assert!(solve_balance_point(&load, &cap, -25.0, 18.0, 0.0).is_err());
        // Non-finite endpoint.
        assert!(solve_balance_point(&load, &cap, f64::NAN, 18.0, 1e-9).is_err());
    }

    #[test]
    fn invalid_line_parameters_are_rejected() {
        assert!(LoadLine::new(18.0, 0.0).is_err());
        assert!(LoadLine::new(18.0, -1.0).is_err());
        assert!(LoadLine::new(f64::INFINITY, 0.5).is_err());
        assert!(CapacityLine::new(8.3, 0.0, 0.25).is_err());
        assert!(CapacityLine::new(8.3, 10.0, 0.0).is_err());
        assert!(CapacityLine::new(8.3, 10.0, -0.25).is_err());
    }

    #[test]
    fn serde_round_trips() {
        let load = LoadLine::new(18.0, 0.5).unwrap();
        let cap = CapacityLine::new(8.3, 10.0, 0.25).unwrap();
        let bp = solve_balance_point(&load, &cap, -25.0, 18.0, 1e-12).unwrap();
        let json = serde_json::to_string(&bp).unwrap();
        let back: BalancePoint = serde_json::from_str(&json).unwrap();
        assert_eq!(bp, back);
    }
}
