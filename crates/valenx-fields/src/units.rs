//! SI dimension-aware units.
//!
//! Every numeric quantity in a `Field`, `ScalarRecord`, or result
//! manifest carries a [`Units`] value. Arithmetic on `Units` checks
//! that dimensions combine correctly — multiplying a velocity by a
//! time yields a length; multiplying two pressures yields an
//! impossible `Pa²` which the code allows (algebra is symmetric) but
//! which display code surfaces honestly.
//!
//! The type follows RFC 0004. Powers of the seven SI base
//! dimensions are stored in the order **L, M, T, I, Θ, N, J** —
//! length, mass, time, electric current, temperature, amount of
//! substance, luminous intensity.
//!
//! Most common units are provided as `pub const` so consumers can
//! write `units::PASCAL` / `units::METER_PER_SECOND` without
//! constructing them by hand.

use std::fmt;
use std::ops::{Div, Mul};

use serde::{Deserialize, Serialize};

/// A physical unit expressed as powers of the 7 SI base dimensions,
/// plus a display hint.
///
/// `display` is a non-essential human-readable symbol. It's skipped
/// during (de)serialization because:
/// - `&'static str` conflicts with serde's `Deserialize<'de>` lifetime
/// - `display_string()` can always synthesise a symbol from `dim`
/// - Canonical constants below supply the symbol at runtime
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Units {
    /// Signed integer powers of `[L, M, T, I, Θ, N, J]`.
    pub dim: [i8; 7],
    /// Linear scale relative to the base SI unit with the same
    /// dimensions (e.g. km → `1000.0`, millisecond → `0.001`).
    pub scale: f64,
    /// Short human-readable symbol (`"Pa"`, `"m/s"`). Optional; when
    /// `None`, renderers fall back to [`Units::display_string`] which
    /// synthesises a symbol from `dim`.
    #[serde(skip, default)]
    pub display: Option<&'static str>,
}

impl Units {
    /// Construct a units value. Most callers should use the constants
    /// at the module bottom and combine them with `*` / `/`.
    pub const fn new(dim: [i8; 7], scale: f64, display: Option<&'static str>) -> Self {
        Self {
            dim,
            scale,
            display,
        }
    }

    /// True if `self` and `other` describe the same physical quantity
    /// (same dimensions) — scale may differ.
    pub fn is_compatible(&self, other: &Units) -> bool {
        self.dim == other.dim
    }

    /// Scale factor needed to convert a value in `self` to a value in
    /// `to`, or `None` when the dimensions disagree.
    ///
    /// ```
    /// # use valenx_fields::units::{Units, METER, KILOMETER};
    /// assert_eq!(KILOMETER.convert_factor(&METER), Some(1000.0));
    /// ```
    pub fn convert_factor(&self, to: &Units) -> Option<f64> {
        if self.is_compatible(to) {
            Some(self.scale / to.scale)
        } else {
            None
        }
    }

    /// Human-readable symbol, synthesised from `dim` when no explicit
    /// `display` is set.
    pub fn display_string(&self) -> String {
        if let Some(s) = self.display {
            return s.to_string();
        }
        const NAMES: [&str; 7] = ["m", "kg", "s", "A", "K", "mol", "cd"];
        let mut num = String::new();
        let mut den = String::new();
        for (i, &p) in self.dim.iter().enumerate() {
            match p.cmp(&0) {
                std::cmp::Ordering::Greater => {
                    push_term(&mut num, NAMES[i], p);
                }
                std::cmp::Ordering::Less => {
                    push_term(&mut den, NAMES[i], -p);
                }
                std::cmp::Ordering::Equal => {}
            }
        }
        match (num.is_empty(), den.is_empty()) {
            (true, true) => "1".to_string(),
            (false, true) => num,
            (true, false) => format!("1/{den}"),
            (false, false) => format!("{num}/{den}"),
        }
    }
}

fn push_term(dst: &mut String, name: &str, power: i8) {
    if !dst.is_empty() {
        dst.push('·');
    }
    dst.push_str(name);
    if power != 1 {
        use std::fmt::Write;
        let _ = write!(dst, "^{power}");
    }
}

impl fmt::Display for Units {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.display_string())
    }
}

impl PartialEq for Units {
    fn eq(&self, other: &Self) -> bool {
        self.dim == other.dim && self.scale == other.scale
    }
}

impl Mul for Units {
    type Output = Units;
    fn mul(self, rhs: Units) -> Units {
        let mut dim = [0i8; 7];
        for (out, (a, b)) in dim.iter_mut().zip(self.dim.iter().zip(rhs.dim.iter())) {
            *out = a.saturating_add(*b);
        }
        Units::new(dim, self.scale * rhs.scale, None)
    }
}

impl Div for Units {
    type Output = Units;
    fn div(self, rhs: Units) -> Units {
        let mut dim = [0i8; 7];
        for (out, (a, b)) in dim.iter_mut().zip(self.dim.iter().zip(rhs.dim.iter())) {
            *out = a.saturating_sub(*b);
        }
        Units::new(dim, self.scale / rhs.scale, None)
    }
}

// ---------------------------------------------------------------------------
// Canonical constants
// ---------------------------------------------------------------------------

/// Dimensionless quantity (ratio, coefficient, count).
pub const DIMENSIONLESS: Units = Units::new([0, 0, 0, 0, 0, 0, 0], 1.0, Some(""));

// Base SI units.
/// SI metre — base unit of length.
pub const METER: Units = Units::new([1, 0, 0, 0, 0, 0, 0], 1.0, Some("m"));
/// SI kilogram — base unit of mass.
pub const KILOGRAM: Units = Units::new([0, 1, 0, 0, 0, 0, 0], 1.0, Some("kg"));
/// SI second — base unit of time.
pub const SECOND: Units = Units::new([0, 0, 1, 0, 0, 0, 0], 1.0, Some("s"));
/// SI ampere — base unit of electric current.
pub const AMPERE: Units = Units::new([0, 0, 0, 1, 0, 0, 0], 1.0, Some("A"));
/// SI kelvin — base unit of thermodynamic temperature.
pub const KELVIN: Units = Units::new([0, 0, 0, 0, 1, 0, 0], 1.0, Some("K"));
/// SI mole — base unit of amount of substance.
pub const MOLE: Units = Units::new([0, 0, 0, 0, 0, 1, 0], 1.0, Some("mol"));
/// SI candela — base unit of luminous intensity.
pub const CANDELA: Units = Units::new([0, 0, 0, 0, 0, 0, 1], 1.0, Some("cd"));

// Common length scales.
/// Kilometre (`1 km = 1000 m`).
pub const KILOMETER: Units = Units::new([1, 0, 0, 0, 0, 0, 0], 1000.0, Some("km"));
/// Centimetre (`1 cm = 0.01 m`).
pub const CENTIMETER: Units = Units::new([1, 0, 0, 0, 0, 0, 0], 0.01, Some("cm"));
/// Millimetre (`1 mm = 0.001 m`).
pub const MILLIMETER: Units = Units::new([1, 0, 0, 0, 0, 0, 0], 0.001, Some("mm"));
/// Micrometre (`1 µm = 1e-6 m`).
pub const MICROMETER: Units = Units::new([1, 0, 0, 0, 0, 0, 0], 1.0e-6, Some("μm"));

// Common time scales.
/// Millisecond (`1 ms = 1e-3 s`).
pub const MILLISECOND: Units = Units::new([0, 0, 1, 0, 0, 0, 0], 1.0e-3, Some("ms"));
/// Microsecond (`1 µs = 1e-6 s`).
pub const MICROSECOND: Units = Units::new([0, 0, 1, 0, 0, 0, 0], 1.0e-6, Some("μs"));
/// Minute (`1 min = 60 s`).
pub const MINUTE: Units = Units::new([0, 0, 1, 0, 0, 0, 0], 60.0, Some("min"));
/// Hour (`1 h = 3600 s`).
pub const HOUR: Units = Units::new([0, 0, 1, 0, 0, 0, 0], 3600.0, Some("h"));

// Common mass scales.
/// Gram (`1 g = 1e-3 kg`).
pub const GRAM: Units = Units::new([0, 1, 0, 0, 0, 0, 0], 1.0e-3, Some("g"));
/// Tonne (`1 t = 1000 kg`).
pub const TONNE: Units = Units::new([0, 1, 0, 0, 0, 0, 0], 1000.0, Some("t"));

// Derived.
/// Velocity in metres per second.
pub const METER_PER_SECOND: Units = Units::new([1, 0, -1, 0, 0, 0, 0], 1.0, Some("m/s"));
/// Force in newtons (`1 N = 1 kg·m/s²`).
pub const NEWTON: Units = Units::new([1, 1, -2, 0, 0, 0, 0], 1.0, Some("N"));
/// Pressure in pascals (`1 Pa = 1 N/m²`).
pub const PASCAL: Units = Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0, Some("Pa"));
/// Pressure in kilopascals (`1 kPa = 1000 Pa`).
pub const KILOPASCAL: Units = Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0e3, Some("kPa"));
/// Pressure in megapascals (`1 MPa = 1e6 Pa`).
pub const MEGAPASCAL: Units = Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0e6, Some("MPa"));
/// Pressure in bar (`1 bar = 1e5 Pa`).
pub const BAR: Units = Units::new([-1, 1, -2, 0, 0, 0, 0], 1.0e5, Some("bar"));
/// Energy in joules (`1 J = 1 N·m`).
pub const JOULE: Units = Units::new([2, 1, -2, 0, 0, 0, 0], 1.0, Some("J"));
/// Power in watts (`1 W = 1 J/s`).
pub const WATT: Units = Units::new([2, 1, -3, 0, 0, 0, 0], 1.0, Some("W"));
/// Frequency in hertz (`1 Hz = 1 s⁻¹`).
pub const HERTZ: Units = Units::new([0, 0, -1, 0, 0, 0, 0], 1.0, Some("Hz"));
/// Electric potential in volts (`1 V = 1 W/A`).
pub const VOLT: Units = Units::new([2, 1, -3, -1, 0, 0, 0], 1.0, Some("V"));
/// Electrical resistance in ohms (`1 Ω = 1 V/A`).
pub const OHM: Units = Units::new([2, 1, -3, -2, 0, 0, 0], 1.0, Some("Ω"));
/// Electric charge in coulombs (`1 C = 1 A·s`).
pub const COULOMB: Units = Units::new([0, 0, 1, 1, 0, 0, 0], 1.0, Some("C"));
/// Capacitance in farads (`1 F = 1 C/V`).
pub const FARAD: Units = Units::new([-2, -1, 4, 2, 0, 0, 0], 1.0, Some("F"));
/// Kinematic viscosity (m²/s).
pub const KINEMATIC_VISCOSITY: Units = Units::new([2, 0, -1, 0, 0, 0, 0], 1.0, Some("m²/s"));
/// Dynamic viscosity (Pa·s).
pub const DYNAMIC_VISCOSITY: Units = Units::new([-1, 1, -1, 0, 0, 0, 0], 1.0, Some("Pa·s"));
/// Mass density (kg/m³).
pub const DENSITY: Units = Units::new([-3, 1, 0, 0, 0, 0, 0], 1.0, Some("kg/m³"));
/// Temperature in degrees Celsius. Shares dimensions with [`KELVIN`];
/// note that the offset conversion (`°C = K − 273.15`) is *not* a
/// linear scale, so use a dedicated helper when converting values.
pub const CELSIUS: Units = Units::new([0, 0, 0, 0, 1, 0, 0], 1.0, Some("°C"));
// Note: Celsius and Kelvin share dimensions. Offset conversion (+273.15)
// is not a linear scale; use a dedicated helper when converting values.

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compatibility() {
        assert!(METER.is_compatible(&KILOMETER));
        assert!(PASCAL.is_compatible(&BAR));
        assert!(!METER.is_compatible(&SECOND));
        assert!(!PASCAL.is_compatible(&NEWTON));
    }

    #[test]
    fn convert_factor_between_scales() {
        assert_eq!(KILOMETER.convert_factor(&METER), Some(1000.0));
        assert_eq!(METER.convert_factor(&KILOMETER), Some(0.001));
        assert_eq!(MINUTE.convert_factor(&SECOND), Some(60.0));
        assert!(PASCAL.convert_factor(&SECOND).is_none());
    }

    #[test]
    fn multiplication_adds_dimensions() {
        // velocity × time = length
        let length = METER_PER_SECOND * SECOND;
        assert!(length.is_compatible(&METER));
    }

    #[test]
    fn division_subtracts_dimensions() {
        // length / time = velocity
        let v = METER / SECOND;
        assert!(v.is_compatible(&METER_PER_SECOND));
    }

    #[test]
    fn derived_units_round_trip() {
        // force / area = pressure
        let area = METER * METER;
        let pressure = NEWTON / area;
        assert!(pressure.is_compatible(&PASCAL));
    }

    #[test]
    fn display_string_synthesises_when_no_symbol() {
        let custom = METER * KELVIN;
        // No explicit display on the product; should synthesise.
        let s = custom.display_string();
        assert!(s.contains('m') && s.contains('K'), "got {s}");
    }

    #[test]
    fn display_string_uses_explicit_symbol() {
        assert_eq!(PASCAL.display_string(), "Pa");
        assert_eq!(METER_PER_SECOND.display_string(), "m/s");
    }

    #[test]
    fn dimensionless_reports_correctly() {
        let ratio = METER / METER;
        assert!(ratio.is_compatible(&DIMENSIONLESS));
    }
}
