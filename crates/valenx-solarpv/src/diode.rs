//! The single-diode photovoltaic cell model.
//!
//! ## The equation
//!
//! A solar cell is modelled as a light-driven current source `Iph` in
//! parallel with a diode, a shunt resistance `Rsh`, and in series with a
//! resistance `Rs`. Kirchhoff's current law at the terminal gives the
//! canonical *five-parameter* single-diode equation
//!
//! ```text
//! I = Iph - I0 * (exp(q * (V + I*Rs) / (n * k * T)) - 1) - (V + I*Rs) / Rsh
//! ```
//!
//! where the *thermal voltage* `Vt = k*T/q` (see
//! [`crate::constants::thermal_voltage`]) collapses the exponent to
//! `(V + I*Rs) / (n * Vt)`.
//!
//! ## Ideal vs. general
//!
//! When `Rs = 0` and `Rsh = +inf` the equation is *explicit* in `I`:
//!
//! ```text
//! I(V) = Iph - I0 * (exp(V / (n * Vt)) - 1)
//! ```
//!
//! and the open-circuit voltage and short-circuit current have closed
//! forms (see [`SingleDiode::voc`] and [`SingleDiode::isc`]). When
//! `Rs > 0` the `I*Rs` term makes the equation *implicit*, so the
//! terminal current at a fixed voltage is found by a damped Newton
//! iteration in [`SingleDiode::current_at`].

use crate::constants::thermal_voltage;
use crate::error::{Result, SolarPvError};

/// A single-diode photovoltaic cell (or series-connected module treated
/// as one lumped cell) at one operating temperature.
///
/// Construct via [`SingleDiode::new`] (validates every parameter) or
/// [`SingleDiode::ideal`] (the explicit `Rs = 0`, `Rsh = inf` case).
/// All currents are in amperes, voltages in volts, resistances in ohms,
/// temperature in kelvin.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SingleDiode {
    /// Photo-generated (light) current `Iph`, in amperes. Proportional to
    /// irradiance; at short circuit in the ideal model `Isc = Iph`.
    pub iph: f64,
    /// Diode reverse-saturation (dark) current `I0`, in amperes. Sets the
    /// recombination floor; typically `1e-12` to `1e-9` A for silicon.
    pub i0: f64,
    /// Diode ideality factor `n` (dimensionless), nominally in `[1, 2]`
    /// for a single junction. Larger `n` softens the knee of the curve.
    pub n: f64,
    /// Absolute cell temperature `T`, in kelvin. Must be strictly
    /// positive.
    pub temperature_k: f64,
    /// Series resistance `Rs`, in ohms (`>= 0`). Models contact and bulk
    /// resistance; `0` recovers the explicit ideal model.
    pub rs: f64,
    /// Shunt resistance `Rsh`, in ohms (`> 0`, may be `+inf`). Models
    /// leakage across the junction; `+inf` recovers the ideal model.
    pub rsh: f64,
}

impl SingleDiode {
    /// Build a fully general single-diode cell, validating every
    /// parameter.
    ///
    /// # Errors
    ///
    /// Returns [`SolarPvError::Invalid`] if any value is non-finite or
    /// out of its physical domain: `iph < 0`, `i0 < 0`, `n` outside
    /// `(0, 10]`, `temperature_k <= 0`, `rs < 0`, or `rsh <= 0`.
    /// (`rsh = +inf` is accepted and denotes "no shunt leakage".)
    pub fn new(iph: f64, i0: f64, n: f64, temperature_k: f64, rs: f64, rsh: f64) -> Result<Self> {
        if !iph.is_finite() || iph < 0.0 {
            return Err(SolarPvError::invalid(
                "iph",
                format!("photocurrent must be finite and >= 0, got {iph}"),
            ));
        }
        if !i0.is_finite() || i0 < 0.0 {
            return Err(SolarPvError::invalid(
                "i0",
                format!("saturation current must be finite and >= 0, got {i0}"),
            ));
        }
        if !n.is_finite() || n <= 0.0 || n > 10.0 {
            return Err(SolarPvError::invalid(
                "n",
                format!("ideality factor must be finite and in (0, 10], got {n}"),
            ));
        }
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(SolarPvError::invalid(
                "temperature_k",
                format!("absolute temperature must be finite and > 0, got {temperature_k}"),
            ));
        }
        if !rs.is_finite() || rs < 0.0 {
            return Err(SolarPvError::invalid(
                "rs",
                format!("series resistance must be finite and >= 0, got {rs}"),
            ));
        }
        // Rsh may legitimately be +inf (ideal, no shunt). Reject only
        // NaN and non-positive values.
        if rsh.is_nan() || rsh <= 0.0 {
            return Err(SolarPvError::invalid(
                "rsh",
                format!("shunt resistance must be > 0 (or +inf), got {rsh}"),
            ));
        }
        Ok(Self {
            iph,
            i0,
            n,
            temperature_k,
            rs,
            rsh,
        })
    }

    /// Build an *ideal* single-diode cell: `Rs = 0`, `Rsh = +inf`.
    ///
    /// In this regime [`current_at`](Self::current_at) is the explicit
    /// `Iph - I0*(exp(V/(n*Vt)) - 1)`, and [`voc`](Self::voc) /
    /// [`isc`](Self::isc) have closed forms.
    ///
    /// # Errors
    ///
    /// Same domain checks as [`new`](Self::new) for `iph`, `i0`, `n`,
    /// `temperature_k`.
    pub fn ideal(iph: f64, i0: f64, n: f64, temperature_k: f64) -> Result<Self> {
        Self::new(iph, i0, n, temperature_k, 0.0, f64::INFINITY)
    }

    /// Thermal voltage `Vt = k*T/q` at this cell's temperature, in volts.
    ///
    /// # Errors
    ///
    /// Cannot fail for a constructed cell (temperature was validated
    /// `> 0`), but propagates [`SolarPvError::Invalid`] for symmetry with
    /// [`crate::constants::thermal_voltage`].
    pub fn thermal_voltage(&self) -> Result<f64> {
        thermal_voltage(self.temperature_k)
    }

    /// `true` when this cell is the explicit ideal case
    /// (`Rs == 0` and `Rsh == +inf`).
    pub fn is_ideal(&self) -> bool {
        self.rs == 0.0 && self.rsh.is_infinite()
    }

    /// Terminal current `I` (amperes) at terminal voltage `voltage_v`
    /// (volts).
    ///
    /// For an ideal cell this evaluates the explicit equation directly.
    /// For `Rs > 0` (or finite `Rsh`) the implicit equation is solved by
    /// a damped Newton iteration whose derivative is known in closed
    /// form, so convergence is quadratic near the root.
    ///
    /// The returned current may be negative (the cell is being driven
    /// past `Voc` into the diode-forward quadrant); callers interested
    /// only in the power-producing quadrant clamp at the relevant
    /// bounds themselves.
    ///
    /// # Errors
    ///
    /// Returns [`SolarPvError::Invalid`] if `voltage_v` is non-finite, or
    /// [`SolarPvError::NoConvergence`] if the Newton solve (general case)
    /// fails to reach tolerance within the iteration budget.
    pub fn current_at(&self, voltage_v: f64) -> Result<f64> {
        if !voltage_v.is_finite() {
            return Err(SolarPvError::invalid(
                "voltage_v",
                format!("voltage must be finite, got {voltage_v}"),
            ));
        }
        let vt = self.thermal_voltage()?;
        let n_vt = self.n * vt;

        // Explicit ideal path: no Rs, no shunt term.
        if self.is_ideal() {
            return Ok(self.iph - self.i0 * ((voltage_v / n_vt).exp() - 1.0));
        }

        // General implicit path. Solve
        //   f(I) = Iph - I0*(exp((V + I*Rs)/(n*Vt)) - 1) - (V + I*Rs)/Rsh - I = 0
        // by Newton's method. f'(I) = -I0*Rs/(n*Vt)*exp(...) - Rs/Rsh - 1,
        // which is strictly negative, so the iteration is well-behaved.
        let shunt_conductance = if self.rsh.is_infinite() {
            0.0
        } else {
            1.0 / self.rsh
        };

        // Seed with the ideal-model current at this voltage (ignoring the
        // resistive drop) — a good starting point near the operating
        // range.
        let mut i = self.iph - self.i0 * ((voltage_v / n_vt).exp() - 1.0);
        if !i.is_finite() {
            i = self.iph;
        }

        const MAX_ITERS: u32 = 200;
        const TOL: f64 = 1e-12;
        let mut residual = f64::INFINITY;
        for _ in 0..MAX_ITERS {
            let v_drop = voltage_v + i * self.rs;
            let exp_term = (v_drop / n_vt).exp();
            let f = self.iph - self.i0 * (exp_term - 1.0) - v_drop * shunt_conductance - i;
            let df = -self.i0 * self.rs / n_vt * exp_term - self.rs * shunt_conductance - 1.0;
            // df is bounded away from zero (<= -1), so this is safe.
            let step = f / df;
            i -= step;
            residual = f.abs();
            if residual <= TOL {
                return Ok(i);
            }
        }
        Err(SolarPvError::no_convergence(
            "current_at",
            MAX_ITERS,
            residual,
        ))
    }

    /// Instantaneous power `P = V * I` (watts) at terminal voltage
    /// `voltage_v` (volts).
    ///
    /// # Errors
    ///
    /// Propagates any error from [`current_at`](Self::current_at).
    pub fn power_at(&self, voltage_v: f64) -> Result<f64> {
        Ok(voltage_v * self.current_at(voltage_v)?)
    }

    /// Open-circuit voltage `Voc` (volts) — the voltage at which the
    /// terminal current is zero.
    ///
    /// For the ideal cell this is the closed form
    /// `Voc = n * Vt * ln(Iph / I0 + 1)`. For the general case it is
    /// found by a bracketed bisection on `current_at(V) = 0`.
    ///
    /// # Errors
    ///
    /// Returns [`SolarPvError::Invalid`] if `i0 == 0` (no dark current →
    /// `Voc` diverges) or `iph == 0` (no light → `Voc = 0`, returned
    /// directly rather than as an error). Returns
    /// [`SolarPvError::NoConvergence`] if the general-case bracket search
    /// fails.
    pub fn voc(&self) -> Result<f64> {
        if self.iph == 0.0 {
            // No illumination: the cell sits at the origin, Voc = 0.
            return Ok(0.0);
        }
        if self.i0 == 0.0 {
            return Err(SolarPvError::invalid(
                "i0",
                "open-circuit voltage diverges when the saturation current is zero",
            ));
        }
        let vt = self.thermal_voltage()?;
        let n_vt = self.n * vt;

        if self.is_ideal() {
            // I(V)=0 => exp(V/n_vt) - 1 = Iph/I0 => V = n_vt*ln(Iph/I0 + 1).
            return Ok(n_vt * (self.iph / self.i0 + 1.0).ln());
        }

        // General case: bracket and bisect current_at(V) == 0. The ideal
        // Voc is an upper bound (series/shunt only ever lower it), and 0
        // gives a positive current (≈ Iph), so [0, voc_ideal] brackets
        // the root.
        let voc_ideal = n_vt * (self.iph / self.i0 + 1.0).ln();
        self.solve_zero_current(0.0, voc_ideal)
    }

    /// Short-circuit current `Isc` (amperes) — the terminal current when
    /// the voltage is zero.
    ///
    /// For the ideal cell this is exactly `Iph`. For the general case it
    /// is `current_at(0.0)`, which with `Rs > 0` is implicit and very
    /// slightly below `Iph`.
    ///
    /// # Errors
    ///
    /// Propagates any error from [`current_at`](Self::current_at).
    pub fn isc(&self) -> Result<f64> {
        self.current_at(0.0)
    }

    /// Bisection helper: find the voltage in `[lo, hi]` where
    /// `current_at(V) == 0`, assuming `current_at(lo) >= 0 >=
    /// current_at(hi)`.
    fn solve_zero_current(&self, lo: f64, hi: f64) -> Result<f64> {
        let mut a = lo;
        let mut b = hi;
        let mut fa = self.current_at(a)?;
        const MAX_ITERS: u32 = 200;
        const TOL: f64 = 1e-12;
        let mut mid = 0.5 * (a + b);
        for _ in 0..MAX_ITERS {
            mid = 0.5 * (a + b);
            let fm = self.current_at(mid)?;
            if fm.abs() <= TOL || (b - a) <= TOL {
                return Ok(mid);
            }
            // Root is where current crosses from + to -.
            if (fa > 0.0) == (fm > 0.0) {
                a = mid;
                fa = fm;
            } else {
                b = mid;
            }
        }
        Err(SolarPvError::no_convergence(
            "voc",
            MAX_ITERS,
            self.current_at(mid)?.abs(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::STC_TEMPERATURE_K;

    /// A representative crystalline-silicon cell at STC.
    fn si_cell() -> SingleDiode {
        // Iph ≈ 3.8 A, I0 ≈ 1e-9 A, n = 1.2, with small Rs and large Rsh.
        SingleDiode::new(3.8, 1.0e-9, 1.2, STC_TEMPERATURE_K, 0.01, 200.0).unwrap()
    }

    fn ideal_cell() -> SingleDiode {
        SingleDiode::ideal(3.8, 1.0e-9, 1.2, STC_TEMPERATURE_K).unwrap()
    }

    #[test]
    fn constructor_rejects_bad_parameters() {
        assert!(SingleDiode::new(-1.0, 1e-9, 1.0, 300.0, 0.0, 100.0).is_err());
        assert!(SingleDiode::new(3.0, -1e-9, 1.0, 300.0, 0.0, 100.0).is_err());
        assert!(SingleDiode::new(3.0, 1e-9, 0.0, 300.0, 0.0, 100.0).is_err());
        assert!(SingleDiode::new(3.0, 1e-9, 11.0, 300.0, 0.0, 100.0).is_err());
        assert!(SingleDiode::new(3.0, 1e-9, 1.0, 0.0, 0.0, 100.0).is_err());
        assert!(SingleDiode::new(3.0, 1e-9, 1.0, 300.0, -0.1, 100.0).is_err());
        assert!(SingleDiode::new(3.0, 1e-9, 1.0, 300.0, 0.0, 0.0).is_err());
        assert!(SingleDiode::new(3.0, 1e-9, 1.0, 300.0, 0.0, f64::NAN).is_err());
        // +inf shunt is allowed.
        assert!(SingleDiode::new(3.0, 1e-9, 1.0, 300.0, 0.0, f64::INFINITY).is_ok());
    }

    /// VALIDATE: I(0) ~ Isc, and for the ideal cell I(0) == Iph exactly.
    #[test]
    fn ideal_short_circuit_current_equals_iph() {
        let c = ideal_cell();
        let i0v = c.current_at(0.0).unwrap();
        assert!((i0v - c.iph).abs() < 1e-12, "I(0) = {i0v}");
        let isc = c.isc().unwrap();
        assert!((isc - c.iph).abs() < 1e-12, "Isc = {isc}");
    }

    /// VALIDATE: I(0) ~ Isc for the general cell — within a hair of Iph
    /// because the resistive drop at V=0 is only I*Rs with small Rs.
    #[test]
    fn general_short_circuit_current_close_to_iph() {
        let c = si_cell();
        let isc = c.isc().unwrap();
        // Rs is tiny, so Isc is just below Iph but very close.
        assert!(isc <= c.iph + 1e-9, "Isc = {isc}");
        assert!((isc - c.iph).abs() < 1e-2, "Isc = {isc}, Iph = {}", c.iph);
    }

    /// VALIDATE: I(Voc) == 0 for both the ideal and general cells.
    #[test]
    fn current_is_zero_at_voc() {
        for c in [ideal_cell(), si_cell()] {
            let voc = c.voc().unwrap();
            let i_at_voc = c.current_at(voc).unwrap();
            assert!(i_at_voc.abs() < 1e-9, "I(Voc) = {i_at_voc} for {c:?}");
        }
    }

    /// Ground truth: the ideal Voc matches the closed form
    /// n*Vt*ln(Iph/I0 + 1) recomputed independently.
    #[test]
    fn ideal_voc_matches_closed_form() {
        let c = ideal_cell();
        let vt = c.thermal_voltage().unwrap();
        let expected = c.n * vt * (c.iph / c.i0 + 1.0).ln();
        let voc = c.voc().unwrap();
        assert!((voc - expected).abs() < 1e-12, "Voc = {voc}");
        // Sanity: a silicon cell's Voc lands near ~0.6-0.65 V.
        assert!((0.55..0.75).contains(&voc), "Voc = {voc}");
    }

    /// The general-case Voc is never above the ideal Voc (series/shunt
    /// resistance can only reduce it).
    #[test]
    fn general_voc_not_above_ideal_voc() {
        let ideal = ideal_cell();
        let real = si_cell();
        let voc_ideal = ideal.voc().unwrap();
        let voc_real = real.voc().unwrap();
        assert!(
            voc_real <= voc_ideal + 1e-9,
            "voc_real = {voc_real}, voc_ideal = {voc_ideal}"
        );
    }

    /// Cross-check the general Newton solver against the explicit ideal
    /// formula: with Rs=0, Rsh=inf the implicit path must reproduce the
    /// closed form to machine precision. Force the implicit path by using
    /// a near-zero (but nonzero) Rs and huge Rsh and comparing.
    #[test]
    fn newton_path_agrees_with_explicit_at_small_rs() {
        let explicit = ideal_cell();
        let near_ideal = SingleDiode::new(3.8, 1.0e-9, 1.2, STC_TEMPERATURE_K, 1e-9, 1e12).unwrap();
        for v in [0.0, 0.1, 0.3, 0.5, 0.55] {
            let ie = explicit.current_at(v).unwrap();
            let ig = near_ideal.current_at(v).unwrap();
            assert!((ie - ig).abs() < 1e-4, "v={v}: ie={ie}, ig={ig}");
        }
    }

    /// The I-V curve is monotonically decreasing in V over [0, Voc].
    #[test]
    fn iv_curve_is_monotonic_decreasing() {
        let c = si_cell();
        let voc = c.voc().unwrap();
        let mut prev = f64::INFINITY;
        for k in 0..=50 {
            let v = voc * (k as f64) / 50.0;
            let i = c.current_at(v).unwrap();
            assert!(i <= prev + 1e-9, "non-monotone at v={v}: {i} > {prev}");
            prev = i;
        }
    }

    /// Power is zero at both endpoints (V=0 → P=0; V=Voc → I=0 → P=0) and
    /// strictly positive somewhere in between.
    #[test]
    fn power_vanishes_at_endpoints_and_is_positive_inside() {
        let c = si_cell();
        let voc = c.voc().unwrap();
        assert!(c.power_at(0.0).unwrap().abs() < 1e-12);
        assert!(c.power_at(voc).unwrap().abs() < 1e-9);
        let p_mid = c.power_at(0.5 * voc).unwrap();
        assert!(p_mid > 0.0, "P(Voc/2) = {p_mid}");
    }

    #[test]
    fn current_at_rejects_non_finite_voltage() {
        let c = si_cell();
        assert!(c.current_at(f64::NAN).is_err());
        assert!(c.current_at(f64::INFINITY).is_err());
    }

    #[test]
    fn voc_zero_when_no_illumination() {
        let dark = SingleDiode::ideal(0.0, 1e-9, 1.0, 300.0).unwrap();
        assert_eq!(dark.voc().unwrap(), 0.0);
    }

    #[test]
    fn voc_errors_when_saturation_current_zero() {
        let c = SingleDiode::ideal(3.0, 0.0, 1.0, 300.0).unwrap();
        assert!(c.voc().is_err());
    }
}
