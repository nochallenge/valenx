//! Temperature-unit helpers.
//!
//! All thermistor physics in this crate works in **kelvin** because the
//! Steinhart-Hart and beta equations are written in terms of absolute
//! temperature. These free functions convert to and from degrees
//! Celsius for callers who think in everyday units. They are pure and
//! total (no validation): the Celsius<->Kelvin map is affine and
//! defined for every finite input.

/// The offset between the Celsius and Kelvin scales, in kelvin.
///
/// `0 degrees C` corresponds to `273.15 K` by definition of the
/// Celsius scale.
pub const KELVIN_AT_ZERO_CELSIUS: f64 = 273.15;

/// Convert a temperature in degrees Celsius to kelvin.
///
/// Adds [`KELVIN_AT_ZERO_CELSIUS`]; the scales share the same degree
/// size, so this is a pure offset.
///
/// ```
/// use valenx_thermistor::units::celsius_to_kelvin;
/// assert!((celsius_to_kelvin(25.0) - 298.15).abs() < 1e-9);
/// ```
pub fn celsius_to_kelvin(celsius: f64) -> f64 {
    celsius + KELVIN_AT_ZERO_CELSIUS
}

/// Convert a temperature in kelvin to degrees Celsius.
///
/// Subtracts [`KELVIN_AT_ZERO_CELSIUS`]; the exact inverse of
/// [`celsius_to_kelvin`].
///
/// ```
/// use valenx_thermistor::units::kelvin_to_celsius;
/// assert!((kelvin_to_celsius(310.15) - 37.0).abs() < 1e-9);
/// ```
pub fn kelvin_to_celsius(kelvin: f64) -> f64 {
    kelvin - KELVIN_AT_ZERO_CELSIUS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_celsius_is_freezing_point() {
        assert!((celsius_to_kelvin(0.0) - 273.15).abs() < 1e-12);
    }

    #[test]
    fn round_trips_celsius() {
        for c in [-40.0, 0.0, 25.0, 37.0, 100.0, 125.0] {
            let back = kelvin_to_celsius(celsius_to_kelvin(c));
            assert!((back - c).abs() < 1e-9, "round trip failed for {c}");
        }
    }

    #[test]
    fn body_temperature() {
        // 37 C is the canonical human body temperature.
        assert!((celsius_to_kelvin(37.0) - 310.15).abs() < 1e-9);
    }
}
